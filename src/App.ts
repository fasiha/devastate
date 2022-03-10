import * as Plot from "@observablehq/plot";
import {configureStore, createSlice} from '@reduxjs/toolkit'
import {createElement as ce, useEffect, useState} from 'react';
import {useDispatch, useSelector} from 'react-redux';

import {get, Question, QuestionBlock} from './questions';

type Score = {
  results: {[k: string]: {result: boolean|undefined, confidence: number|undefined}},
  show: boolean
};
const initialScore: Score = {
  results: {},
  show: false
};
const CONFIDENCES = [55, 65, 75, 85, 95];
// These are my results from when I read the book :D
const me: [boolean, number][] = [
  [true, 95], [false, 55], [true, 65],  [true, 85], [true, 65], [true, 85], [false, 55], [true, 65],
  [true, 55], [false, 85], [true, 95],  [true, 75], [true, 55], [true, 95], [true, 85],  [true, 95],
  [true, 95], [true, 55],  [false, 55], [true, 75], [true, 55], [true, 85], [true, 65],  [true, 65],
  [true, 85], [true, 55],  [true, 55],  [true, 55], [true, 95], [true, 55], [true, 95],  [true, 75],
  [true, 65], [true, 55],  [true, 95],  [true, 85], [true, 85], [true, 95], [true, 95],  [true, 75]
];

const slice = createSlice({
  name: 'score',
  initialState: initialScore,
  reducers: {
    append: (state, action) => ({...state, results: {...state.results, ...action.payload}}),
    done: state => ({...state, show: true}),
    fill: state => ({
      ...state,
      results: Object.fromEntries(
          Object.keys(state.results).map((key, i) => [key, {result: me[i][0], confidence: me[i][1]}]))
    }),
  }
});
export const store = configureStore({reducer: {score: slice.reducer}});
type RootState = ReturnType<typeof store.getState>;
const actions = slice.actions;
const totalAnswered = (s: RootState) => ({
  total: Object.keys(s.score.results).length,
  answered: Object.values(s.score.results).filter(x => x.result !== undefined && x.confidence !== undefined).length
});
const show = (s: RootState) => s.score.show;
const almostDone = (s: RootState) => {
  const {total, answered} = totalAnswered(s);
  return total - answered <= 3;
};

function hash(question: Question): string { return (question.question || '') + question.options.join(''); }

interface SummaryProps {
  detailed: boolean;
}
function Summary({detailed}: SummaryProps) {
  const {answered, total} = useSelector(totalAnswered);
  const showState = useSelector(show);
  const resultState = useSelector((s: RootState) => s.score.results);

  const dispatch = useDispatch();

  const plotData: {x: number, y: number}[] = [];

  useEffect(() => {
    if (detailed) {
      const dot = Plot.dot(plotData, {x: 'x', y: 'y', stroke: 'red', r: 10});
      const link = Plot.link([1], {x1: 50, y1: 50, x2: 100, y2: 100, strokeOpacity: 0.2});
      document.querySelector('#plot')?.replaceChildren(Plot.plot({
        marks: [dot, link],
        grid: true,
        x: {label: 'When you feel ░░% sure of your answer…', ticks: CONFIDENCES},
        y: {label: '… you\'re right ░░% of the time '},
        style: {background: "black", color: "white"},
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

    plotData.push({x: conf, y: successes / total * 100});
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
  const dispatch = useDispatch();

  const choice: number|undefined = resultState.result === undefined ? undefined
                                   : resultState.result             ? question.answer
                                                                    : Number(!question.answer);

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
               dispatch(actions.append({[unique]: {...resultState, result: num === question.answer}}));
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
  const dispatch = useDispatch();
  useEffect(() => {
    get().then(list => {
      setQuestionBlocks(list);
      dispatch(actions.append(Object.fromEntries(list.flatMap(l => l.questions.map(q => [hash(q), {}])))))
    });
  }, []);
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