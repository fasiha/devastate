import {configureStore, createSlice} from '@reduxjs/toolkit'
import React, {createElement as ce, useEffect, useState} from 'react';
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

const slice = createSlice({
  name: 'score',
  initialState: initialScore,
  reducers: {
    append: (state, action) => ({...state, results: {...state.results, ...action.payload}}),
    done: state => ({...state, show: true}),
    fill: state => ({
      ...state,
      results:
          Object.fromEntries(Object.keys(state.results)
                                 .map(key => [key, {result: Math.random() < 0.5, confidence: pickrand(CONFIDENCES)}]))
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

function hash(question: Question): string { return (question.question || '') + question.options.join(''); }
function pickrand<T>(v: T[]): T { return v[Math.floor(Math.random() * v.length)]; }

function Summary() {
  const {answered, total} = useSelector(totalAnswered);
  const dispatch = useDispatch();
  if (answered !== total) {
    return ce('p', {}, `${answered} of ${total} question(s) answered!`, ' ',
              ce('button', {onClick: () => dispatch(actions.fill())}, 'Random?'));
  }
  return ce('p', {}, ce('button', {onClick: () => dispatch(actions.done())}, 'Show results?!'));
}

interface QProps {
  question: Question;
}
function Q({question}: QProps) {
  // `unique` is a "hash" of the question, so even if somehow we rerender, it'll be the same
  const unique = hash(question);
  const selector = useSelector<RootState, Score['results'][string]>(s => s.score.results[unique]) || {};
  const dispatch = useDispatch();

  const choice: number|undefined = selector.result === undefined ? undefined
                                   : selector.result             ? question.answer
                                                                 : Number(!question.answer);

  function makeAnswers(content: string, num: number) {
    const radioGroup = unique + 'question';
    // <input type="radio" id="dewey" name="drone" value="dewey">
    // <label for="dewey">Dewey</label>
    const id = radioGroup + num;
    return [
      ce('input', {
        type: 'radio',
        id,
        name: radioGroup,
        checked: choice === num,
        onChange: e => {
          if (e.target.value) { dispatch(actions.append({[unique]: {...selector, result: num === question.answer}})); }
        }
      }),
      ce('label', {htmlFor: id}, content)
    ];
  };
  const pairs = question.options.flatMap((option, num) => makeAnswers(option, num));

  function makeConfidences(content: string, num: number) {
    const radioGroup = unique + 'conf';
    const id = radioGroup + num;
    const thisConf = CONFIDENCES[num];
    return [
      ce('input', {
        type: 'radio',
        id,
        name: radioGroup,
        checked: selector.confidence === thisConf,
        onChange: e => {
          if (e.target.value) { dispatch(actions.append({[unique]: {...selector, confidence: thisConf}})); }
        }
      }),
      ce('label', {htmlFor: id}, content)
    ];
  };
  const confidences = CONFIDENCES.flatMap((option, num) => makeConfidences(`${option}%`, num))

  return ce('li', {}, (question.question || '') + ' ', ...pairs, ' — confidence: ', ...confidences);
}

interface BlockProps {
  block: QuestionBlock;
}
function Block({block}: BlockProps) {
  return ce('div', {}, block.title, ce('ul', {}, ...block.questions.map(question => ce(Q, {question}))))
}

function App() {
  const [questionBlocks, setQuestionBlocks] = useState<QuestionBlock[]|undefined>(undefined);
  const dispatch = useDispatch();
  useEffect(() => {
    get().then(list => {
      setQuestionBlocks(list);
      for (const b of list) {
        for (const q of b.questions) {
          dispatch(actions.append({[hash(q)]: {result: undefined, confidence: undefined}}));
        }
      }
    });
  }, []);
  if (questionBlocks) { return ce('div', null, ce(Summary), ...questionBlocks.map(block => ce(Block, {block}))); }
  return ce('div', null, 'Waiting for data!');
}

export default App;