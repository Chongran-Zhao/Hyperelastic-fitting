import React, { useState } from 'react';

const ModelViewer = () => {
  const [modelType, setModelType] = useState('gaussian');
  const [parameters, setParameters] = useState({
    gaussian: { mu: 1.0 },
    nonGaussian: { mu: 1.0, N: 10 }
  });

  const handleModelChange = (e) => {
    setModelType(e.target.value);
  };

  const handleParameterChange = (param, value) => {
    setParameters(prev => ({
      ...prev,
      [modelType]: {
        ...prev[modelType],
        [param]: parseFloat(value)
      }
    }));
  };

  return (
    <div className="model-viewer">
      <div className="model-selector">
        <h2>微观力学模型参数配置</h2>
        <select value={modelType} onChange={handleModelChange}>
          <option value="gaussian">Gaussian 模型</option>
          <option value="nonGaussian">Non-Gaussian 模型</option>
        </select>
      </div>

      <div className="parameters-panel">
        <h3>模型参数</h3>
        {modelType === 'gaussian' ? (
          <div className="parameter-group">
            <label>
              μ (剪切模量):
              <input
                type="number"
                value={parameters.gaussian.mu}
                onChange={(e) => handleParameterChange('mu', e.target.value)}
                step="0.1"
              />
            </label>
          </div>
        ) : (
          <div className="parameter-group">
            <label>
              μ (剪切模量):
              <input
                type="number"
                value={parameters.nonGaussian.mu}
                onChange={(e) => handleParameterChange('mu', e.target.value)}
                step="0.1"
              />
            </label>
            <label>
              N (链段数):
              <input
                type="number"
                value={parameters.nonGaussian.N}
                onChange={(e) => handleParameterChange('N', e.target.value)}
                step="1"
                min="1"
              />
            </label>
          </div>
        )}
      </div>

      <div className="results-panel">
        <h3>计算结果</h3>
        {/* 这里将添加计算结果的可视化组件 */}
      </div>
    </div>
  );
};

export default ModelViewer;