import React, {createElement as ce} from "../_snowpack/pkg/react.js";
import ReactDOM from "../_snowpack/pkg/react-dom.js";
import {Provider} from "../_snowpack/pkg/react-redux.js";
import App, {store} from "./App.js";
ReactDOM.render(ce(React.StrictMode, null, ce(Provider, {store}, ce(App))), document.getElementById("root"));
