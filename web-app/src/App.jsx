import { useState } from 'react'
import './App.css'
import ModelViewer from './components/ModelViewer'

function App() {
  return (
    <div className="app">
      <header className="app-header">
        <h1>超弹性材料微观力学模型</h1>
      </header>
      <main className="app-main">
        <ModelViewer />
      </main>
    </div>
  )
}

export default App
