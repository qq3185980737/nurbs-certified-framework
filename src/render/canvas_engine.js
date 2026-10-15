// ============================================================================
// Canvas 渲染引擎 —— UI 使用的规范渲染器（由 index.html 内联副本迁出，单一事实来源）
// 由 build.mjs 打包注入 index.html；改动此处后运行 `node build.mjs` 同步到浏览器。
// 依赖：math_core(求值/AABB) + solver(GeometryState)
// ============================================================================
// 依赖：wasm_backend(求值/AABB/GeometryState，转发到 C++ WASM 数学核)
import { evaluatePoint, computeAABB, GeometryState } from '../wasm_backend.js';

const COLORS = {
    background: '#1a1a2e', grid: '#2d2d44', gridMajor: '#3d3d5c', axis: '#4a4a6a',
    curve1: '#00d9ff', curve2: '#ff6b6b', controlPolygon1: 'rgba(0, 217, 255, 0.3)',
    controlPolygon2: 'rgba(255, 107, 107, 0.3)', controlPoint1: '#00a8cc', controlPoint2: '#ff5252',
    intersection: '#ffd700', text: '#e0e0e0'
};

class Viewport {
    constructor(width, height) {
        this.width = width; this.height = height;
        this.zoom = 1.0;
        this.offsetX = width / 2; this.offsetY = height / 2;
    }
    
    worldToScreen(wx, wy) {
        return { x: wx * this.zoom + this.offsetX, y: -wy * this.zoom + this.offsetY };
    }
    
    screenToWorld(sx, sy) {
        return { x: (sx - this.offsetX) / this.zoom, y: -(sy - this.offsetY) / this.zoom };
    }
    
    fitToBounds(bounds, padding = 0.1) {
        const rangeX = bounds.maxX - bounds.minX || 1;
        const rangeY = bounds.maxY - bounds.minY || 1;
        const scaleX = (this.width * (1 - 2 * padding)) / rangeX;
        const scaleY = (this.height * (1 - 2 * padding)) / rangeY;
        this.zoom = Math.min(scaleX, scaleY);
        const centerX = (bounds.minX + bounds.maxX) / 2;
        const centerY = (bounds.minY + bounds.maxY) / 2;
        this.offsetX = this.width / 2 - centerX * this.zoom;
        this.offsetY = this.height / 2 + centerY * this.zoom;
    }
}

class CanvasEngine {
    constructor(canvas) {
        this.canvas = canvas;
        this.ctx = canvas.getContext('2d');
        this.viewport = new Viewport(canvas.width, canvas.height);
        this._rafId = null;
        this.setupEventListeners();
        this._applyHiDPI();
    }
    
    // HiDPI：canvas 后台缓冲 = CSS 像素 × devicePixelRatio，绘制坐标保持 CSS 像素
    _applyHiDPI() {
        const dpr = (typeof window !== 'undefined' && window.devicePixelRatio) || 1;
        if (dpr <= 1) return;
        const rect = this.canvas.getBoundingClientRect();
        if (!rect.width || !rect.height) return;
        const w = Math.round(rect.width * dpr);
        const h = Math.round(rect.height * dpr);
        if (this.canvas.width !== w) this.canvas.width = w;
        if (this.canvas.height !== h) this.canvas.height = h;
        this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    }
    
    // rAF 合并重绘：滚轮/拖拽等高频事件在一帧内只重绘一次，避免每次事件全量渲染卡顿
    requestRender() {
        if (this._rafId != null) return;
        this._rafId = requestAnimationFrame(() => {
            this._rafId = null;
            this.render();
        });
    }
    
    setupEventListeners() {
        const canvas = this.canvas;
        let isDragging = false, lastX = 0, lastY = 0;
        
        canvas.addEventListener('wheel', (e) => {
            e.preventDefault();
            const delta = e.deltaY > 0 ? 0.9 : 1.1;
            this.viewport.zoom = Math.max(0.01, Math.min(100, this.viewport.zoom * delta));
            this.requestRender();
        }, { passive: false });
        
        canvas.addEventListener('mousedown', (e) => { isDragging = true; lastX = e.clientX; lastY = e.clientY; });
        canvas.addEventListener('mousemove', (e) => {
            if (isDragging) {
                this.viewport.offsetX += e.clientX - lastX;
                this.viewport.offsetY += e.clientY - lastY;
                lastX = e.clientX; lastY = e.clientY;
                this.requestRender();
            }
        });
        canvas.addEventListener('mouseup', () => { isDragging = false; });
        
        window.addEventListener('resize', () => {
            const rect = canvas.getBoundingClientRect();
            this._applyHiDPI();
            this.viewport.width = rect.width; this.viewport.height = rect.height;
            this.viewport.offsetX = rect.width / 2; this.viewport.offsetY = rect.height / 2;
            this.requestRender();
        });
    }
    
    clear() {
        this.ctx.fillStyle = COLORS.background;
        this.ctx.fillRect(0, 0, this.viewport.width, this.viewport.height);
    }
    
    drawGrid() {
        const ctx = this.ctx, vp = this.viewport;
        ctx.strokeStyle = COLORS.grid; ctx.lineWidth = 1;
        const topLeft = vp.screenToWorld(0, 0);
        const bottomRight = vp.screenToWorld(vp.width, vp.height);
        const gridSize = vp.zoom > 5 ? 0.5 : vp.zoom > 2 ? 1 : vp.zoom < 0.5 ? 5 : 1;
        
        ctx.beginPath();
        for (let x = Math.floor(topLeft.x / gridSize) * gridSize; x <= bottomRight.x; x += gridSize) {
            const s = vp.worldToScreen(x, 0);
            ctx.moveTo(s.x, 0); ctx.lineTo(s.x, vp.height);
        }
        for (let y = Math.floor(bottomRight.y / gridSize) * gridSize; y <= topLeft.y; y += gridSize) {
            const s = vp.worldToScreen(0, y);
            ctx.moveTo(0, s.y); ctx.lineTo(vp.width, s.y);
        }
        ctx.stroke();
        
        ctx.strokeStyle = COLORS.axis; ctx.lineWidth = 2; ctx.beginPath();
        const xAxis = vp.worldToScreen(0, 0);
        ctx.moveTo(0, xAxis.y); ctx.lineTo(vp.width, xAxis.y);
        const yAxis = vp.worldToScreen(0, 0);
        ctx.moveTo(yAxis.x, 0); ctx.lineTo(yAxis.x, vp.height);
        ctx.stroke();
    }
    
    drawCurve(curve, color) {
        const ctx = this.ctx, vp = this.viewport;
        const points = [];
        const steps = Math.max(50, curve.controlPoints.length * 20);
        for (let i = 0; i <= steps; i++) {
            const u = curve.uMin + (curve.uMax - curve.uMin) * i / steps;
            try { points.push(evaluatePoint(curve, u)); } catch(e) {}
        }
        if (points.length === 0) return;
        
        ctx.beginPath();
        const start = vp.worldToScreen(points[0].x, points[0].y);
        ctx.moveTo(start.x, start.y);
        for (let i = 1; i < points.length; i++) {
            const s = vp.worldToScreen(points[i].x, points[i].y);
            ctx.lineTo(s.x, s.y);
        }
        ctx.strokeStyle = color; ctx.lineWidth = 2.5; ctx.lineCap = 'round'; ctx.lineJoin = 'round'; ctx.stroke();
    }
    
    drawControlPolygon(curve, color) {
        const ctx = this.ctx, vp = this.viewport;
        ctx.beginPath();
        const start = vp.worldToScreen(curve.controlPoints[0].x, curve.controlPoints[0].y);
        ctx.moveTo(start.x, start.y);
        for (let i = 1; i < curve.controlPoints.length; i++) {
            const s = vp.worldToScreen(curve.controlPoints[i].x, curve.controlPoints[i].y);
            ctx.lineTo(s.x, s.y);
        }
        ctx.strokeStyle = color; ctx.lineWidth = 1; ctx.setLineDash([5, 5]); ctx.stroke(); ctx.setLineDash([]);
    }
    
    drawControlPoints(curve, color) {
        const ctx = this.ctx, vp = this.viewport;
        for (let i = 0; i < curve.controlPoints.length; i++) {
            const pt = curve.controlPoints[i];
            const s = vp.worldToScreen(pt.x, pt.y);
            ctx.beginPath(); ctx.arc(s.x, s.y, 6, 0, Math.PI * 2);
            ctx.fillStyle = color; ctx.fill();
            ctx.strokeStyle = '#fff'; ctx.lineWidth = 2; ctx.stroke();
        }
    }
    
    drawResult(result) {
        const ctx = this.ctx, vp = this.viewport;
        const s1 = vp.worldToScreen(result.point1.x, result.point1.y);
        
        ctx.beginPath(); ctx.arc(s1.x, s1.y, 10, 0, Math.PI * 2);
        ctx.fillStyle = 'rgba(255, 215, 0, 0.3)'; ctx.fill();
        ctx.beginPath(); ctx.arc(s1.x, s1.y, 6, 0, Math.PI * 2);
        ctx.fillStyle = COLORS.intersection; ctx.fill();
        ctx.strokeStyle = '#fff'; ctx.lineWidth = 2; ctx.stroke();
        
        if (result.state === GeometryState.SEPARATED && result.point2) {
            const s2 = vp.worldToScreen(result.point2.x, result.point2.y);
            ctx.beginPath(); ctx.moveTo(s1.x, s1.y); ctx.lineTo(s2.x, s2.y);
            ctx.strokeStyle = '#a55eea'; ctx.lineWidth = 2; ctx.setLineDash([5, 3]); ctx.stroke(); ctx.setLineDash([]);
        }
    }
    
    render(curve1 = null, curve2 = null, results = []) {
        this.clear();
        this.drawGrid();
        if (curve1) {
            this.drawControlPolygon(curve1, COLORS.controlPolygon1);
            this.drawCurve(curve1, COLORS.curve1);
            this.drawControlPoints(curve1, COLORS.controlPoint1);
        }
        if (curve2) {
            this.drawControlPolygon(curve2, COLORS.controlPolygon2);
            this.drawCurve(curve2, COLORS.curve2);
            this.drawControlPoints(curve2, COLORS.controlPoint2);
        }
        for (const result of results) { this.drawResult(result); }
    }
    
    fitCurves(curve1, curve2) {
        const box1 = computeAABB(curve1), box2 = computeAABB(curve2);
        const combinedBox = {
            minX: Math.min(box1.minX, box2.minX),
            minY: Math.min(box1.minY, box2.minY),
            maxX: Math.max(box1.maxX, box2.maxX),
            maxY: Math.max(box1.maxY, box2.maxY)
        };
        this.viewport.fitToBounds(combinedBox, 0.15);
        this.render(curve1, curve2, []);
    }
}

// ============================================================================
// 导出
// ============================================================================
export { CanvasEngine, Viewport, COLORS };
