import * as Plot from "../_snowpack/pkg/@observablehq/plot.js";
import {configureStore, createSlice} from "../_snowpack/pkg/@reduxjs/toolkit.js";
import {createElement as ce, useEffect, useState} from "../_snowpack/pkg/react.js";
import {useDispatch, useSelector} from "../_snowpack/pkg/react-redux.js";
import {lookups} from "./confintervals.js";
import {get} from "./questions.js";
const initialScore = {
  results: [],
  show: false
};
const CONFIDENCES = [55, 65, 75, 85, 95];
const me = [
  [true, 95, 1],
  [false, 55, 1],
  [true, 65, 1],
  [true, 85, 0],
  [true, 65, 1],
  [true, 85, 1],
  [false, 55, 1],
  [true, 65, 0],
  [true, 55, 1],
  [false, 85, 1],
  [true, 95, 1],
  [true, 75, 1],
  [true, 55, 0],
  [true, 95, 0],
  [true, 85, 1],
  [true, 95, 1],
  [true, 95, 0],
  [true, 55, 1],
  [false, 55, 0],
  [true, 75, 0],
  [true, 55, 0],
  [true, 85, 0],
  [true, 65, 0],
  [true, 65, 0],
  [true, 85, 0],
  [true, 55, 0],
  [true, 55, 1],
  [true, 55, 0],
  [true, 95, 0],
  [true, 55, 1],
  [true, 95, 1],
  [true, 75, 0],
  [true, 65, 1],
  [true, 55, 0],
  [true, 95, 1],
  [true, 85, 1],
  [true, 85, 0],
  [true, 95, 1],
  [true, 95, 0],
  [true, 75, 0]
];
const confidenceToIdx = new Map(CONFIDENCES.map((c, i) => [c, i]));
function serializeResults(results) {
  let suggestion = "v1-" + results.map((o) => `${o.choice ?? "x"}${confidenceToIdx.get(o.confidence ?? -1) ?? "x"}`).join("");
  suggestion = suggestion.replace(/x+$/, "");
  return suggestion;
}
function deserializeResults(s, qs) {
  const ret = [];
  let pieces = [];
  if (s.startsWith("v1-")) {
    pieces = s.slice(3).match(/(.{2}|.{1})/g) || [];
  }
  for (const [i, q] of qs.entries()) {
    const piece = pieces[i];
    if (piece) {
      const [choice, confIdx] = piece.split("").map((o) => parseInt(o));
      ret.push({
        choice: isNaN(choice) ? void 0 : choice,
        confidence: isNaN(confIdx) ? void 0 : CONFIDENCES[confIdx],
        result: isNaN(choice) ? void 0 : choice === q.answer
      });
    } else {
      ret.push({choice: void 0, confidence: void 0, result: void 0});
    }
  }
  return ret;
}
function lerp(x1, x2, y1, y2, x) {
  const mu = (x - x1) / (x2 - x1);
  return y1 * (1 - mu) + y2 * mu;
}
function pureArrayReplaceIdx(v, newElement, idx) {
  const ret = v.slice();
  ret[idx] = newElement;
  return ret;
}
const slice = createSlice({
  name: "score",
  initialState: initialScore,
  reducers: {
    append: function(state, action) {
      return {...state, results: pureArrayReplaceIdx(state.results, action.payload[1], action.payload[0])};
    },
    set: function(state, action) {
      return {...state, results: action.payload};
    },
    done: (state) => ({...state, show: true}),
    fill: (state) => ({...state, results: me.map((m) => ({result: m[0], confidence: m[1], choice: m[2]}))})
  }
});
export const store = configureStore({reducer: {score: slice.reducer}});
const useAppDispatch = () => useDispatch();
const actions = slice.actions;
const totalAnswered = (s) => ({
  total: Object.keys(s.score.results).length,
  answered: Object.values(s.score.results).filter((x) => x.result !== void 0 && x.confidence !== void 0).length
});
const show = (s) => s.score.show;
const almostDone = (s) => {
  const {total, answered} = totalAnswered(s);
  return total - answered <= 3;
};
const serializedHash = (s) => serializeResults(s.score.results);
function Summary({detailed}) {
  const {answered, total} = useSelector(totalAnswered);
  const showState = useSelector(show);
  const resultState = useSelector((s) => s.score.results);
  const dispatch = useAppDispatch();
  const plotData = [];
  useEffect(() => {
    if (detailed) {
      plotData.sort((a, b) => a.x - b.x);
      const cis = [];
      for (const {total: n, x: p} of plotData) {
        for (const [q, ci] of lookups(n, p)) {
          if (q < 40 || q > 60) {
            cis.push({x: p, y: ci / n * 100, q});
          }
        }
      }
      const ciDots = Plot.line(cis, {
        x: "x",
        y: "y",
        z: "q",
        stroke: "q",
        curve: "basis",
        strokeWidth: (d) => lerp(90, 50, 1, 2, 2 * Math.abs(50 - d.q)),
        strokeOpacity: 0.5
      });
      const ciText = Plot.text(cis, Plot.selectLast({
        x: "x",
        y: "y",
        z: "q",
        text: (d) => {
          const q = d.q;
          const ci = 2 * Math.abs(50 - q);
          return `${q > 50 ? "↓" : "↑"}${ci}% CI`;
        },
        textAnchor: "start",
        dx: 10
      }));
      const dot = Plot.dot(plotData, {x: "x", y: "y", stroke: "red", fill: "red", r: 10});
      const link = Plot.link([1], {x1: 50, y1: 50, x2: 100, y2: 100, stroke: "orange", strokeOpacity: 0.5});
      document.querySelector("#plot")?.replaceChildren(Plot.plot({
        marks: [dot, ciDots, ciText, link],
        grid: true,
        x: {label: "When you feel ░░% sure of your answer…", ticks: CONFIDENCES},
        y: {label: "… you're right ░░% of the time "},
        style: {background: "black", color: "white"},
        r: {type: "linear", domain: [0, 1], range: [0, 10]},
        color: {type: "linear", domain: [0, 50, 100], range: ["white", "orange", "white"]}
      }));
    } else if (showState) {
      document.querySelector("#plot")?.scrollIntoView({behavior: "smooth"});
    }
  }, [detailed, showState]);
  if (answered !== total) {
    const fill = answered === 0 ? ce("button", {onClick: () => dispatch(actions.fill())}, "Fill?") : "";
    return ce("p", {}, `${answered} of ${total} question(s) answered! `, fill);
  }
  if (!showState) {
    return ce("p", {}, ce("button", {onClick: () => dispatch(actions.done())}, "Show results?!"));
  }
  if (!detailed) {
    return ce("p");
  }
  const buckets = new Map();
  for (const [idx, {result, confidence}] of Object.values(resultState).entries()) {
    const c = confidence ?? -1;
    const r = !!result;
    buckets.set(c, (buckets.get(c) || []).concat({result: r, idx}));
  }
  const bullets = [];
  for (const [conf, results] of buckets) {
    const successes = results.reduce((p, c) => p + Number(c.result), 0);
    const total2 = results.length;
    const pct = (successes / total2 * 100).toFixed(1);
    const rightIdxs = [];
    const wrongIdxs = [];
    for (const {result, idx} of results) {
      (result ? rightIdxs : wrongIdxs).push(idx + 1);
    }
    bullets.push(`${conf}% ➜ ${pct}% = ${successes}/${total2}: [${rightIdxs.join(" ")}] right vs [${wrongIdxs.join(" ")}] wrong`);
    plotData.push({x: conf, y: successes / total2 * 100, successes, total: total2});
  }
  bullets.sort();
  return ce("div", {}, ce("p", {}, "Results! ", ce("button", {
    onClick: () => {
      try {
        navigator.clipboard.writeText(window.location.href);
      } catch {
        alert("Please copy the URL from the address bar!");
      }
    }
  }, "Copy URL")), ce("ul", {}, ...bullets.map((text) => ce("li", {}, text))));
}
function Q({question}) {
  const unique = "q" + question.idx;
  const resultState = useSelector((s) => s.score.results[question.idx - 1]) || {};
  const showState = useSelector(show);
  const dispatch = useAppDispatch();
  const choice = resultState.choice;
  if (showState && choice !== void 0) {
    const options = question.options.join(" or ");
    const summary = `${question.question || ""} ${options}. You said ${question.options[choice]}, ${resultState.confidence}%. `;
    const result = `${resultState.result ? "✅" : "❌"}!`;
    const comment = question.comment ? ` ${question.comment}` : "";
    return ce("li", {className: "question", value: question.idx}, summary, result, comment);
  }
  function makeAnswers(content, num) {
    const radioGroup = unique + "question";
    const id = radioGroup + num;
    return [
      ce("label", {htmlFor: id}, ce("input", {
        type: "radio",
        id,
        name: radioGroup,
        checked: choice === num,
        onChange: (e) => {
          if (e.target.value) {
            dispatch(actions.append([question.idx - 1, {...resultState, choice: num, result: num === question.answer}]));
          }
        }
      }), " " + content),
      " "
    ];
  }
  ;
  const pairs = question.options.flatMap((option, num) => makeAnswers(option, num));
  function makeConfidences(content, num) {
    const radioGroup = unique + "conf";
    const id = radioGroup + num;
    const thisConf = CONFIDENCES[num];
    return [
      ce("label", {htmlFor: id}, ce("input", {
        type: "radio",
        id,
        name: radioGroup,
        checked: resultState.confidence === thisConf,
        onChange: (e) => {
          if (e.target.value) {
            dispatch(actions.append([question.idx - 1, {...resultState, confidence: thisConf}]));
          }
        }
      }), " " + content),
      " "
    ];
  }
  ;
  const confidences = CONFIDENCES.flatMap((option, num) => makeConfidences(`${option}%`, num));
  const className = "question " + (resultState.confidence !== void 0 && resultState.result !== void 0 ? "answered" : "unanswered");
  return ce("li", {className, value: question.idx}, ...question.question ? [question.question, ce("br")] : [""], ...pairs, ce("br"), "Confidence: ", ...confidences);
}
function Block({block}) {
  return ce("div", {}, block.title, ce("ol", {}, ...block.questions.map((question) => ce(Q, {question}))));
}
function App() {
  const [questionBlocks, setQuestionBlocks] = useState(void 0);
  const almost = useSelector(almostDone);
  const dispatch = useAppDispatch();
  useEffect(() => {
    get().then((list) => {
      setQuestionBlocks(list);
      const deserialized = deserializeResults(window.location.hash.slice(1), list.flatMap((o) => o.questions));
      if (deserialized) {
        dispatch(actions.set(deserialized));
      } else {
        const empty = list.flatMap((l) => l.questions.map((q) => ({result: void 0, confidence: void 0, choice: void 0})));
        dispatch(actions.set(empty));
      }
    });
  }, []);
  const serialized = useSelector(serializedHash);
  useEffect(() => {
    if (questionBlocks) {
      window.location.hash = serialized;
    }
  }, [serialized]);
  if (questionBlocks) {
    return ce("div", {className: almost ? "almost-done" : void 0}, ce(Summary, {detailed: true}), ...questionBlocks.map((block) => ce(Block, {block})), ce(Summary, {detailed: false}));
  }
  return ce("div", null, "Waiting for data!");
}
export default App;
