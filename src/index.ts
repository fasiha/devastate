import React, {createElement as ce} from 'react';
import ReactDOM from 'react-dom';

import App from './App.jsx';

ReactDOM.render(ce(React.StrictMode, null, ce(App)), document.getElementById('root'));