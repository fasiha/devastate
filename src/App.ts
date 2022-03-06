import React, {createElement as ce, useEffect, useState} from 'react';

import {get, Question, QuestionBlock} from './questions';

const SCORE: {[k: string]: {result: boolean|undefined, confidence: number|undefined}} = {};
const CONFIDENCES = [55, 65, 75, 85, 95];

interface QProps {
  question: Question;
}
function Q({question}: QProps) {
  const [choice, setChoice] = useState<number|undefined>(undefined);
  const [conf, setConf] = useState<number|undefined>(undefined);

  // `unique` is a "hash" of the question, so even if somehow we rerender, it'll be the same
  const unique = (question.question || '') + question.options.join('');
  // enroll it in the scoreboard immediately
  if (!(unique in SCORE)) { SCORE[unique] = {result: undefined, confidence: undefined}; }

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
          if (e.target.value) {
            setChoice(num);
            SCORE[unique].result = num === question.answer;
            console.log(SCORE[unique])
          }
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
        checked: conf === thisConf,
        onChange: e => {
          if (e.target.value) {
            setConf(thisConf);
            SCORE[unique].confidence = thisConf;
            console.log(SCORE[unique])
          }
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
  useEffect(() => { get().then(list => { setQuestionBlocks(list); }); }, []);
  if (questionBlocks) { return ce('div', null, ...questionBlocks.map(block => ce(Block, {block}))); }
  return ce('div', null, 'Waiting for data!');
}

export default App;