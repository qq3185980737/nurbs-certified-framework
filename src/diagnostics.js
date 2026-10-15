/**
 * Diagnostics Module - 全局错误监控与诊断增强
 * 
 * 提供应用启动时的健康检查、错误捕获和 UI 提示功能。
 */

'use strict';

/**
 * 初始化全局错误监控
 */
export function initDiagnostics() {
    const errorContainer = createErrorContainer();
    
    // 捕获未处理的异常
    window.onerror = function(message, source, lineno, colno, error) {
        const errorInfo = {
            type: 'JavaScript Error',
            message: message,
            source: source,
            line: lineno,
            column: colno,
            stack: error ? error.stack : 'N/A'
        };
        displayError(errorInfo, errorContainer);
        console.error('[DIAGNOSTICS] Global Error:', errorInfo);
        return false;
    };
    
    // 捕获未处理的 Promise 拒绝
    window.onunhandledrejection = function(event) {
        const errorInfo = {
            type: 'Unhandled Promise Rejection',
            reason: event.reason,
            message: event.reason ? event.reason.message || String(event.reason) : 'Unknown reason'
        };
        displayError(errorInfo, errorContainer);
        console.error('[DIAGNOSTICS] Unhandled Rejection:', errorInfo);
    };
    
    console.log('[DIAGNOSTICS] Global error monitoring initialized');
}

/**
 * 创建错误通知容器
 */
function createErrorContainer() {
    let container = document.getElementById('error-notification-container');
    
    if (!container) {
        container = document.createElement('div');
        container.id = 'error-notification-container';
        container.style.cssText = `
            position: fixed;
            top: 0;
            left: 0;
            right: 0;
            z-index: 9999;
            pointer-events: none;
        `;
        document.body.appendChild(container);
    }
    
    return container;
}

/**
 * 在 UI 上显示错误信息
 */
function displayError(errorInfo, container) {
    const alertBox = document.createElement('div');
    alertBox.style.cssText = `
        background: linear-gradient(135deg, #cc5555, #ff6b6b);
        color: #fff;
        padding: 12px 20px;
        margin: 10px;
        border-radius: 6px;
        box-shadow: 0 4px 12px rgba(255, 107, 107, 0.4);
        font-family: 'Segoe UI', sans-serif;
        font-size: 13px;
        pointer-events: auto;
        cursor: pointer;
        animation: slideDown 0.3s ease-out;
    `;
    
    alertBox.innerHTML = `
        <div style="font-weight: 700; margin-bottom: 6px;">⚠️ ${errorInfo.type}</div>
        <div style="font-family: 'Courier New', monospace; font-size: 12px; opacity: 0.9;">
            ${escapeHtml(errorInfo.message)}
        </div>
        ${errorInfo.source ? `<div style="font-size: 11px; opacity: 0.7; margin-top: 4px;">📍 ${errorInfo.source}:${errorInfo.line}:${errorInfo.column}</div>` : ''}
        <div style="margin-top: 8px; font-size: 11px; opacity: 0.6;">点击关闭</div>
    `;
    
    alertBox.addEventListener('click', () => {
        alertBox.remove();
    });
    
    container.appendChild(alertBox);
    
    // 自动移除（5 秒后）
    setTimeout(() => {
        if (alertBox.parentNode) {
            alertBox.remove();
        }
    }, 10000);
}

/**
 * HTML 转义
 */
function escapeHtml(text) {
    const div = document.createElement('div');
    div.textContent = text;
    return div.innerHTML;
}

/**
 * 模块加载健康检查
 */
export async function checkModuleHealth() {
    const healthReport = {
        modules: {},
        allLoaded: true,
        errors: []
    };
    
    const modulesToCheck = [
        { name: 'math_core', path: './nurbs/math_core.js' },
        { name: 'canvas_engine', path: './render/canvas_engine.js' },
        { name: 'unified_solver', path: './solver/unified_solver.js' }
    ];
    
    for (const mod of modulesToCheck) {
        try {
            await import(mod.path);
            healthReport.modules[mod.name] = { status: 'OK', loaded: true };
            console.log(`[HEALTH CHECK] ✅ ${mod.name} loaded successfully`);
        } catch (error) {
            healthReport.modules[mod.name] = { 
                status: 'FAILED', 
                loaded: false, 
                error: error.message 
            };
            healthReport.allLoaded = false;
            healthReport.errors.push({
                module: mod.name,
                message: error.message
            });
            console.error(`[HEALTH CHECK] ❌ ${mod.name} failed to load:`, error.message);
        }
    }
    
    return healthReport;
}

/**
 * Canvas 上下文检查
 */
export function checkCanvasContext(canvas) {
    if (!canvas) {
        return {
            success: false,
            error: 'Canvas 元素未找到'
        };
    }
    
    const ctx = canvas.getContext('2d');
    
    if (!ctx) {
        return {
            success: false,
            error: '无法获取 Canvas 2D 上下文'
        };
    }
    
    return {
        success: true,
        context: ctx,
        width: canvas.width,
        height: canvas.height
    };
}

/**
 * 生成诊断报告
 */
export function generateDiagnosticReport() {
    const report = {
        timestamp: new Date().toISOString(),
        userAgent: navigator.userAgent,
        screenResolution: `${screen.width}x${screen.height}`,
        devicePixelRatio: window.devicePixelRatio,
        webglSupport: checkWebGLSupport(),
        localStorageAvailable: checkLocalStorage(),
        modulesLoaded: true
    };
    
    return report;
}

function checkWebGLSupport() {
    try {
        const canvas = document.createElement('canvas');
        return !!(window.WebGLRenderingContext && 
                  (canvas.getContext('webgl') || canvas.getContext('experimental-webgl')));
    } catch (e) {
        return false;
    }
}

function checkLocalStorage() {
    try {
        const test = '__localStorage_test__';
        localStorage.setItem(test, test);
        localStorage.removeItem(test);
        return true;
    } catch (e) {
        return false;
    }
}