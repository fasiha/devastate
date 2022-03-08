import React, {createElement as ce} from 'react';
import ReactDOM from 'react-dom';
import {Provider} from 'react-redux';

import App, {store} from './App.jsx';

ReactDOM.render(ce(React.StrictMode, null, ce(Provider, {store}, ce(App))), document.getElementById('root'));