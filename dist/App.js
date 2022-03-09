import * as Plot from "../_snowpack/pkg/@observablehq/plot.js";
import {configureStore, createSlice} from "../_snowpack/pkg/@reduxjs/toolkit.js";
import {createElement as ce, useEffect, useState} from "../_snowpack/pkg/react.js";
import {useDispatch, useSelector} from "../_snowpack/pkg/react-redux.js";
import {get} from "./questions.js";
const initialScore = {
  results: {},
  show: false
};
const CONFIDENCES = [55, 65, 75, 85, 95];
const slice = createSlice({
  name: "score",
  initialState: initialScore,
  reducers: {
    append: (state, action) => ({...state, results: {...state.results, ...action.payload}}),
    done: (state) => ({...state, show: true}),
    fill: (state) => ({
      ...state,
      results: Object.fromEntries(Object.keys(state.results).map((key) => [key, {result: Math.random() < 0.5, confidence: pickrand(CONFIDENCES)}]))
    })
  }
});
export const store = configureStore({reducer: {score: slice.reducer}});
const actions = slice.actions;
const totalAnswered = (s) => ({
  total: Object.keys(s.score.results).length,
  answered: Object.values(s.score.results).filter((x) => x.result !== void 0 && x.confidence !== void 0).length
});
const show = (s) => s.score.show;
function hash(question) {
  return (question.question || "") + question.options.join("");
}
function pickrand(v) {
  return v[Math.floor(Math.random() * v.length)];
}
function Summary() {
  const {answered, total} = useSelector(totalAnswered);
  const showState = useSelector(show);
  const resultState = useSelector((s) => s.score.results);
  const dispatch = useDispatch();
  const plotData = [];
  useEffect(() => {
    if (showState) {
      const dot = Plot.dot(plotData, {x: "x", y: "y"});
      const link = Plot.link([1], {x1: 50, y1: 50, x2: 100, y2: 100, strokeOpacity: 0.2});
      document.querySelector("#plot")?.append(Plot.plot({
        marks: [dot, link],
        grid: true,
        x: {label: "When you feel ░░% sure of your answer…"},
        y: {label: "… you're right ░░% of the time "}
      }));
    }
  });
  if (answered !== total) {
    return ce("p", {}, `${answered} of ${total} question(s) answered!`, " ", ce("button", {onClick: () => dispatch(actions.fill())}, "Random?"));
  }
  if (!showState) {
    return ce("p", {}, ce("button", {onClick: () => dispatch(actions.done())}, "Show results?!"));
  }
  const buckets = new Map();
  for (const {result, confidence} of Object.values(resultState)) {
    const c = confidence ?? -1;
    const r = !!result;
    buckets.set(c, (buckets.get(c) || []).concat(r));
  }
  const bullets = [];
  for (const [conf, results] of buckets) {
    const successes = results.reduce((p, c) => p + Number(c), 0);
    const total2 = results.length;
    const pct = (successes / total2 * 100).toFixed(1);
    bullets.push(`${conf}% ➜ ${pct}% = ${successes}/${total2}`);
    plotData.push({x: conf, y: successes / total2 * 100});
  }
  bullets.sort();
  return ce("div", {}, ce("p", {}, "Results!"), ce("ol", {}, ...bullets.map((text) => ce("li", {}, text))));
}
function Q({question}) {
  const unique = hash(question);
  const resultState = useSelector((s) => s.score.results[unique]) || {};
  const showState = useSelector(show);
  const dispatch = useDispatch();
  const choice = resultState.result === void 0 ? void 0 : resultState.result ? question.answer : Number(!question.answer);
  if (showState && choice !== void 0) {
    const options = question.options.join(" or ");
    const summary = `${question.question || ""} ${options}. You said ${question.options[choice]}. `;
    const result = `${resultState.result ? "✅" : "❌"}!`;
    const comment = question.comment ? ` ${question.comment}` : "";
    return ce("li", {}, summary, result, comment);
  }
  function makeAnswers(content, num) {
    const radioGroup = unique + "question";
    const id = radioGroup + num;
    return [
      ce("input", {
        type: "radio",
        id,
        name: radioGroup,
        checked: choice === num,
        onChange: (e) => {
          if (e.target.value) {
            dispatch(actions.append({[unique]: {...resultState, result: num === question.answer}}));
          }
        }
      }),
      ce("label", {htmlFor: id}, content)
    ];
  }
  ;
  const pairs = question.options.flatMap((option, num) => makeAnswers(option, num));
  function makeConfidences(content, num) {
    const radioGroup = unique + "conf";
    const id = radioGroup + num;
    const thisConf = CONFIDENCES[num];
    return [
      ce("input", {
        type: "radio",
        id,
        name: radioGroup,
        checked: resultState.confidence === thisConf,
        onChange: (e) => {
          if (e.target.value) {
            dispatch(actions.append({[unique]: {...resultState, confidence: thisConf}}));
          }
        }
      }),
      ce("label", {htmlFor: id}, content)
    ];
  }
  ;
  const confidences = CONFIDENCES.flatMap((option, num) => makeConfidences(`${option}%`, num));
  return ce("li", {}, (question.question || "") + " ", ...pairs, " — confidence: ", ...confidences);
}
function Block({block}) {
  return ce("div", {}, block.title, ce("ul", {}, ...block.questions.map((question) => ce(Q, {question}))));
}
function App() {
  const [questionBlocks, setQuestionBlocks] = useState(void 0);
  const dispatch = useDispatch();
  useEffect(() => {
    get().then((list) => {
      setQuestionBlocks(list);
      for (const b of list) {
        for (const q of b.questions) {
          dispatch(actions.append({[hash(q)]: {result: void 0, confidence: void 0}}));
        }
      }
    });
  }, []);
  if (questionBlocks) {
    return ce("div", null, ce(Summary), ...questionBlocks.map((block) => ce(Block, {block})));
  }
  return ce("div", null, "Waiting for data!");
}
export default App;
