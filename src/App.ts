import React, {createElement as ce, useEffect, useState} from 'react';

import {get, Question, QuestionBlock} from './questions';

interface QProps {
  question: Question;
}
function Q({question}: QProps) { return ce('li', {}, question.question); }

interface BlockProps {
  block: QuestionBlock;
}
function Block({block}: BlockProps) {
  return ce('div', {}, block.title, ce('ul', {}, ...block.questions.map(question => ce(Q, {question}))))
}

function App() {
  // Create the count state.
  const [questionBlocks, setQuestionBlocks] = useState<QuestionBlock[]|undefined>(undefined);
  // Update the count (+1 every second).
  useEffect(() => {get().then(x => setQuestionBlocks(x))}, []);
  // Return the App component.
  if (questionBlocks) { return ce('div', null, ...questionBlocks.map(block => ce(Block, {block}))); }
  return ce('div', null, 'Waiting for data!');
}

export default App;