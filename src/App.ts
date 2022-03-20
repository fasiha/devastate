import * as Plot from "@observablehq/plot";
import {configureStore, createSlice, PayloadAction} from '@reduxjs/toolkit'
import {createElement as ce, useEffect, useState} from 'react';
import {useDispatch, useSelector} from 'react-redux';

import {lookups} from "./confintervals";
import {get, hash, Question, QuestionBlock} from './questions';

/*
Business logic types and data
*/
type Score = {
  results: {[k: string]: {result: boolean|undefined, confidence: number|undefined, choice: number|undefined}},
  show: boolean
};
const initialScore: Score = {
  results: {},
  show: false
};
const CONFIDENCES = [55, 65, 75, 85, 95];
// These are my results from when I read the book :D
const me: [boolean, number, number][] = [
  [true, 95, 1], [false, 55, 0], [true, 65, 1],  [true, 85, 0], [true, 65, 1],  [true, 85, 1], [false, 55, 0],
  [true, 65, 0], [true, 55, 1],  [false, 85, 0], [true, 95, 1], [true, 75, 1],  [true, 55, 0], [true, 95, 0],
  [true, 85, 1], [true, 95, 1],  [true, 95, 0],  [true, 55, 1], [false, 55, 1], [true, 75, 0], [true, 55, 0],
  [true, 85, 0], [true, 65, 0],  [true, 65, 0],  [true, 85, 0], [true, 55, 0],  [true, 55, 1], [true, 55, 0],
  [true, 95, 0], [true, 55, 1],  [true, 95, 1],  [true, 75, 0], [true, 65, 1],  [true, 55, 0], [true, 95, 1],
  [true, 85, 1], [true, 85, 0],  [true, 95, 1],  [true, 95, 0], [true, 75, 0]
];

const confidenceToIdx = new Map(CONFIDENCES.map((c, i) => [c, i]))
function serializeResults(results: Score['results']): string {
  console.log('running serialize')
  let suggestion =
      'v1-' +
      Object.values(results).map(o => `${o.choice ?? 'x'}${confidenceToIdx.get(o.confidence ?? -1) ?? 'x'}`).join('');
  // map undefined (unanswered questions/confidences) to 'x' above.
  // Then trim excess 'x's
  suggestion = suggestion.replace(/x+$/, '')
  return suggestion;
}
function deserializeResults(s: string, qs: Question[]): Score['results']|undefined {
  if (s.startsWith('v1-')) {
    const ret: Score['results'] = {};
    // throw away leading version indicator and split into two-char tuples
    const pieces = s.slice(3).match(/../g);
    if (!pieces) { return undefined; }
    for (const [i, piece] of pieces.entries()) {
      const question = qs[i];
      if (!question) { break; } // all done. Garbage in serialized input
      const [choice, confIdx] = piece.split('').map(o => parseInt(o));
      ret[hash(question)] = {
        choice: isNaN(choice) ? undefined : choice,
        confidence: isNaN(confIdx) ? undefined : CONFIDENCES[confIdx],
        result: choice === question.answer
      };
    }
    console.log('returning', ret)
    return ret;
  }
  return undefined;
}

function lerp(x1: number, x2: number, y1: number, y2: number, x: number): number {
  const mu = (x - x1) / (x2 - x1);
  return (y1 * (1 - mu) + y2 * mu);
}

/*
Redux! App state.
*/
const slice = createSlice({
  name: 'score',
  initialState: initialScore,
  reducers: {
    append: (state, action: PayloadAction<Score['results']>) =>
        ({...state, results: {...state.results, ...action.payload}}),
    done: state => ({...state, show: true}),
    fill: state => ({
      ...state,
      results: Object.fromEntries(
          Object.keys(state.results).map((key, i) => [key, {result: me[i][0], confidence: me[i][1], choice: me[i][2]}]))
    }),
    clear: state => ({
      ...state,
      results: Object.fromEntries(
          Object.keys(state.results).map(k => [k, {result: undefined, confidence: undefined, choice: undefined}]))
    }),
  }
});
export const store = configureStore({reducer: {score: slice.reducer}});
type RootState = ReturnType<typeof store.getState>;
const useAppDispatch = () => useDispatch<typeof store.dispatch>();
const actions = slice.actions;
// Reusable selectors
const totalAnswered = (s: RootState) => ({
  total: Object.keys(s.score.results).length,
  answered: Object.values(s.score.results).filter(x => x.result !== undefined && x.confidence !== undefined).length
});
const show = (s: RootState) => s.score.show;
const almostDone = (s: RootState) => {
  const {total, answered} = totalAnswered(s);
  return total - answered <= 3;
};
const serializedHash = (s: RootState) => serializeResults(s.score.results);

/*
React components!
*/
interface SummaryProps {
  detailed: boolean;
}
function Summary({detailed}: SummaryProps) {
  const {answered, total} = useSelector(totalAnswered);
  const showState = useSelector(show);
  const resultState = useSelector((s: RootState) => s.score.results);

  const dispatch = useAppDispatch();

  const plotData: {x: number, y: number, successes: number, total: number}[] = [];

  useEffect(() => {
    if (detailed) {
      plotData.sort((a, b) => a.x - b.x);

      const cis: {x: number, y: number, q: number}[] = [];
      for (const {total: n, x: p} of plotData) {
        for (const [q, ci] of lookups(n, p)) {
          if (q < 40 || q > 60) { cis.push({x: p, y: ci / n * 100, q}); }
        }
      }

      const ciDots = Plot.line(cis, {
        x: 'x',
        y: 'y',
        z: 'q',
        stroke: 'q',
        curve: 'basis',
        strokeWidth: (d: typeof cis[0]) => lerp(90, 50, 1, 2, 2 * Math.abs(50 - d.q)),
        strokeOpacity: 0.5,
      });
      const ciText = Plot.text(cis, Plot.selectLast({
        x: "x",
        y: "y",
        z: "q",
        text: (d: typeof cis[0]) => {
          const q = d.q;
          const ci = 2 * Math.abs(50 - q);
          return `${q > 50 ? '↓' : '↑'}${ci}% CI`;
        },
        textAnchor: "start",
        dx: 10
      }));
      const dot = Plot.dot(plotData, {x: 'x', y: 'y', stroke: 'red', fill: 'red', r: 10});
      const link = Plot.link([1], {x1: 50, y1: 50, x2: 100, y2: 100, stroke: 'orange', strokeOpacity: 0.5});
      document.querySelector('#plot')?.replaceChildren(Plot.plot({
        marks: [dot, ciDots, ciText, link],
        grid: true,
        x: {label: 'When you feel ░░% sure of your answer…', ticks: CONFIDENCES},
        y: {label: '… you\'re right ░░% of the time '},
        style: {background: "black", color: "white"},
        r: {type: "linear", domain: [0, 1], range: [0, 10]},
        color: {type: 'linear', domain: [0, 50, 100], range: ["white", "orange", "white"]},
      }));
    } else if (showState) {
      document.querySelector('#plot')?.scrollIntoView({behavior: 'smooth'});
    }
  }, [detailed, showState]);

  if (answered !== total) {
    const fill = answered === 0 ? ce('button', {onClick: () => dispatch(actions.fill())}, 'Fill?') : '';
    return ce('p', {}, `${answered} of ${total} question(s) answered! `, fill);
  }
  if (!showState) { return ce('p', {}, ce('button', {onClick: () => dispatch(actions.done())}, 'Show results?!')); }

  if (!detailed) { return ce('p'); }

  const buckets: Map<number, {result: boolean, idx: number}[]> = new Map();
  for (const [idx, {result, confidence}] of Object.values(resultState).entries()) {
    const c = confidence ?? -1; // TypeScript pacification
    const r = !!result;         // TypeScript pacification
    buckets.set(c, (buckets.get(c) || []).concat({result: r, idx}))
  }
  const bullets: string[] = [];
  for (const [conf, results] of buckets) {
    const successes = results.reduce((p, c) => p + Number(c.result), 0);
    const total = results.length;
    const pct = (successes / total * 100).toFixed(1);
    const rightIdxs: number[] = [];
    const wrongIdxs: number[] = [];
    for (const {result, idx} of results) { (result ? rightIdxs : wrongIdxs).push(idx + 1); }
    bullets.push(
        `${conf}% ➜ ${pct}% = ${successes}/${total}: [${rightIdxs.join(' ')}] right vs [${wrongIdxs.join(' ')}] wrong`);

    plotData.push({x: conf, y: successes / total * 100, successes, total});
  }
  bullets.sort(); // confidence comes first so we can sort these lexicographically
  return ce(
      'div',
      {},
      ce('p', {}, 'Results!'),
      ce('ul', {}, ...bullets.map(text => ce('li', {}, text))),
  );
}

interface QProps {
  question: Question;
}
function Q({question}: QProps) {
  // `unique` is a "hash" of the question, so even if somehow we rerender, it'll be the same
  const unique = hash(question);
  const resultState = useSelector((s: RootState) => s.score.results[unique]) || {};
  const showState = useSelector(show);
  const dispatch = useAppDispatch();

  const choice = resultState.choice;

  if (showState && choice !== undefined) {
    const options = question.options.join(' or ');
    const summary =
        `${question.question || ''} ${options}. You said ${question.options[choice]}, ${resultState.confidence}%. `;
    const result = `${resultState.result ? '✅' : '❌'}!`;
    const comment = question.comment ? ` ${question.comment}` : '';
    return ce('li', {className: 'question', value: question.idx}, summary, result, comment);
  }

  function makeAnswers(content: string, num: number) {
    const radioGroup = unique + 'question';
    // <input type="radio" id="dewey" name="drone" value="dewey">
    // <label for="dewey">Dewey</label>
    const id = radioGroup + num;
    return [
      ce('label', {htmlFor: id}, ce('input', {
           type: 'radio',
           id,
           name: radioGroup,
           checked: choice === num,
           onChange: e => {
             if (e.target.value) {
               dispatch(actions.append({[unique]: {...resultState, choice: num, result: num === question.answer}}));
             }
           }
         }),
         ' ' + content),
      ' '
    ];
  };
  const pairs = question.options.flatMap((option, num) => makeAnswers(option, num));

  function makeConfidences(content: string, num: number) {
    const radioGroup = unique + 'conf';
    const id = radioGroup + num;
    const thisConf = CONFIDENCES[num];
    return [
      ce('label', {htmlFor: id}, ce('input', {
           type: 'radio',
           id,
           name: radioGroup,
           checked: resultState.confidence === thisConf,
           onChange: e => {
             if (e.target.value) { dispatch(actions.append({[unique]: {...resultState, confidence: thisConf}})); }
           }
         }),
         ' ' + content),
      ' '
    ];
  };
  const confidences = CONFIDENCES.flatMap((option, num) => makeConfidences(`${option}%`, num))

  const className =
      'question ' +
      (resultState.confidence !== undefined && resultState.result !== undefined ? 'answered' : 'unanswered');
  return ce('li', {className, value: question.idx}, ...(question.question ? [question.question, ce('br')] : ['']),
            ...pairs, ce('br'), 'Confidence: ', ...confidences);
}

interface BlockProps {
  block: QuestionBlock;
}
function Block({block}: BlockProps) {
  return ce('div', {}, block.title, ce('ol', {}, ...block.questions.map(question => ce(Q, {question}))))
}

function App() {
  const [questionBlocks, setQuestionBlocks] = useState<QuestionBlock[]|undefined>(undefined);
  const almost = useSelector(almostDone);
  const dispatch = useAppDispatch();
  useEffect(() => {
    get().then(list => {
      setQuestionBlocks(list);
      const empty = Object.fromEntries(list.flatMap(
          l => l.questions.map(q => [hash(q), {result: undefined, confidence: undefined, choice: undefined}])));
      const deserialized = deserializeResults(window.location.hash.slice(1), list.flatMap(o => o.questions)) || {};
      dispatch(actions.append({...empty, ...deserialized}))
    });
  }, []);
  const serialized = useSelector(serializedHash);
  useEffect(() => {
    // don't update the hash too soon!
    if (questionBlocks) { window.location.hash = serialized; }
  }, [serialized]);

  if (questionBlocks) {
    return ce(
        'div',
        {className: almost ? 'almost-done' : undefined},
        ce(Summary, {detailed: true}),
        ...questionBlocks.map(block => ce(Block, {block})),
        ce(Summary, {detailed: false}),
    );
  }
  return ce('div', null, 'Waiting for data!');
}

export default App;