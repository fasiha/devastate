import React, {createElement as ce, useEffect, useState} from 'react';

import {get, Question, QuestionBlock} from './questions';

interface QProps {
  question: Question;
}
function Q({question}: QProps) {
  let q = '';
  if (question.question) { q = question.question + ' '; }
  q += question.options.join(' or ') + '?';
  return ce('li', {}, q);
}

interface BlockProps {
  block: QuestionBlock;
}
function Block({block}: BlockProps) {
  return ce('div', {}, block.title, ce('ul', {}, ...block.questions.map(question => ce(Q, {question}))))
}

function App() {
  const [questionBlocks, setQuestionBlocks] = useState<QuestionBlock[]|undefined>(undefined);
  useEffect(() => { get().then(x => setQuestionBlocks(x)); }, []);
  if (questionBlocks) { return ce('div', null, ...questionBlocks.map(block => ce(Block, {block}))); }
  return ce('div', null, 'Waiting for data!');
}

export default App;