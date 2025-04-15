import React from 'react';
import { HashRouter as Router, Routes, Route } from 'react-router-dom';
import './App.css';
import DataInputPage from './pages/DataInputPage/DataInputPage';
import ResultsPage from './pages/ResultsPage/ResultsPage';

function App() {
  return (
    <Router>
      <div className="App">
        <Routes>
          <Route path="/" element={<DataInputPage />} />
          <Route path="/analysis/:jobId" element={<ResultsPage />} />
        </Routes>
      </div>
    </Router>
  );
}

export default App;
