
将用户提供的文献，根据以下幻灯片模板，制作为组会汇报用幻灯片，使用学术风格的中文编排内容，内容由浅至深，听众为同课题组研究生和导师，因此你的内容不能过于深奥，但也不能过于科普。

需要确保内容足够详细，将文献工作讲解到位，包括其背景、方法、创新点、意义等。保持使用 html 文件呈现内容。


<!DOCTYPE html>
<html lang="zh-CN">

<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>通用学术汇报幻灯片模板</title>
    
    <!-- MathJax 配置 - 支持数学公式 -->
    <script>
        MathJax = {
            tex: {
                inlineMath: [['$', '$'], ['\\(', '\\)']],
                displayMath: [['$$', '$$'], ['\\[', '\\]']]
            },
            svg: {
                fontCache: 'global'
            }
        };
    </script>
    <script type="text/javascript" id="MathJax-script" async
        src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js">
    </script>
    
    <style>
        /* ==================== 主题变量 ==================== */
        :root {
            /* 默认主题 - 专业蓝 */
            --primary-bg: #f8f9fa;
            --secondary-bg: #ffffff;
            --text-color: #2c3e50;
            --header-color: #2563eb;
            --accent-color: #3b82f6;
            --button-bg: #e5e7eb;
            --button-hover-bg: #d1d5db;
            --border-color: #d1d5db;
            --note-bg: #eff6ff;
            --note-border: #3b82f6;
            --shadow-color: rgba(0, 0, 0, 0.1);
            --code-bg: #f3f4f6;
            --code-color: #1f2937;
            --highlight-bg: #fef3c7;
            --success-color: #10b981;
            --warning-color: #f59e0b;
            --error-color: #ef4444;
        }

        /* 深色主题 */
        [data-theme="dark"] {
            --primary-bg: #1a1a1a;
            --secondary-bg: #262626;
            --text-color: #e5e5e5;
            --header-color: #60a5fa;
            --accent-color: #3b82f6;
            --button-bg: #404040;
            --button-hover-bg: #525252;
            --border-color: #404040;
            --note-bg: #1e3a5f;
            --note-border: #60a5fa;
            --shadow-color: rgba(0, 0, 0, 0.5);
            --code-bg: #374151;
            --code-color: #f3f4f6;
            --highlight-bg: #854d0e;
        }

        /* 绿色主题 */
        [data-theme="green"] {
            --primary-bg: #f0fdf4;
            --secondary-bg: #ffffff;
            --text-color: #14532d;
            --header-color: #16a34a;
            --accent-color: #22c55e;
            --button-bg: #dcfce7;
            --button-hover-bg: #bbf7d0;
            --border-color: #bbf7d0;
            --note-bg: #f0fdf4;
            --note-border: #22c55e;
        }

        /* ==================== 全局样式 ==================== */
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        html, body {
            width: 100%;
            height: 100%;
            background: var(--primary-bg);
            color: var(--text-color);
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', 'Roboto', 
                       'Helvetica Neue', 'Microsoft YaHei', 'PingFang SC', sans-serif;
            overflow: hidden;
            transition: background-color 0.3s ease, color 0.3s ease;
        }

        /* ==================== 演示容器 ==================== */
        #presentation-container {
            position: relative;
            width: 100vw;
            height: 56.25vw; /* 16:9 宽高比 */
            max-height: 100vh;
            max-width: 177.78vh;
            margin: auto;
            top: 50%;
            transform: translateY(-50%);
            background: var(--secondary-bg);
            box-shadow: 0 10px 40px var(--shadow-color);
            border-radius: 8px;
            overflow: hidden;
            transition: all 0.3s ease;
        }

        /* ==================== 幻灯片基础样式 ==================== */
        .slide {
            width: 100%;
            height: 100%;
            padding: 3vw;
            display: none;
            flex-direction: column;
            opacity: 0;
            transition: opacity 0.5s ease-in-out;
            position: absolute;
            overflow: hidden;
        }

        .slide.active {
            display: flex;
            opacity: 1;
            z-index: 1;
        }

        .slide-content {
            width: 100%;
            height: 100%;
            overflow-y: auto;
            overflow-x: hidden;
            padding-right: 1vw;
        }

        /* 自定义滚动条 */
        .slide-content::-webkit-scrollbar {
            width: 8px;
        }

        .slide-content::-webkit-scrollbar-track {
            background: var(--button-bg);
            border-radius: 4px;
        }

        .slide-content::-webkit-scrollbar-thumb {
            background: var(--border-color);
            border-radius: 4px;
        }

        .slide-content::-webkit-scrollbar-thumb:hover {
            background: var(--button-hover-bg);
        }

        /* ==================== 标题样式 ==================== */
        h1, h2, h3, h4 {
            color: var(--header-color);
            margin-bottom: 1.5vw;
            font-weight: 600;
            line-height: 1.3;
        }

        h1 {
            font-size: 3.5vw;
            border-bottom: 3px solid var(--accent-color);
            padding-bottom: 1vw;
        }

        h2 {
            font-size: 2.8vw;
            border-bottom: 2px solid var(--accent-color);
            padding-bottom: 0.8vw;
        }

        h3 {
            font-size: 2.2vw;
            color: var(--text-color);
        }

        h4 {
            font-size: 1.8vw;
            color: var(--text-color);
            opacity: 0.9;
        }

        /* ==================== 文本样式 ==================== */
        p {
            font-size: 1.6vw;
            line-height: 1.8;
            margin-bottom: 1.2vw;
        }

        ul, ol {
            margin: 1vw 0;
            padding-left: 2vw;
        }

        li {
            font-size: 1.5vw;
            line-height: 1.8;
            margin-bottom: 0.5vw;
        }

        strong {
            color: var(--header-color);
            font-weight: 600;
        }

        em {
            font-style: italic;
            opacity: 0.9;
        }

        /* ==================== 特殊幻灯片布局 ==================== */
        
        /* 标题页 */
        .title-slide {
            justify-content: center;
            align-items: center;
            text-align: center;
            background: linear-gradient(135deg, var(--secondary-bg) 0%, var(--primary-bg) 100%);
        }

        .title-slide h1 {
            font-size: 4vw;
            border: none;
            margin-bottom: 2vw;
            background: linear-gradient(90deg, var(--header-color), var(--accent-color));
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
        }

        .title-slide .subtitle {
            font-size: 2vw;
            color: var(--text-color);
            opacity: 0.8;
            margin-bottom: 3vw;
        }

        .title-slide .author {
            font-size: 1.8vw;
            margin-bottom: 0.5vw;
        }

        .title-slide .institution {
            font-size: 1.4vw;
            opacity: 0.7;
            margin-bottom: 0.3vw;
        }

        .title-slide .date {
            font-size: 1.3vw;
            opacity: 0.6;
        }

        /* 两列布局 */
        .columns {
            display: flex;
            gap: 2vw;
            width: 100%;
            align-items: flex-start;
        }

        .column {
            flex: 1;
        }

        .column-40 {
            flex: 0 0 40%;
        }

        .column-60 {
            flex: 0 0 60%;
        }

        /* 三列布局 */
        .three-columns {
            display: flex;
            gap: 1.5vw;
            width: 100%;
        }

        .three-columns .column {
            flex: 1;
        }

        /* ==================== 组件样式 ==================== */
        
        /* 图片和图表 */
        figure {
            margin: 1.5vw 0;
            text-align: center;
        }

        img {
            max-width: 100%;
            max-height: 45vh;
            height: auto;
            border-radius: 8px;
            box-shadow: 0 4px 12px var(--shadow-color);
            transition: transform 0.3s ease;
        }

        img:hover {
            transform: scale(1.02);
        }

        figcaption {
            margin-top: 0.8vw;
            font-size: 1.2vw;
            color: var(--text-color);
            opacity: 0.7;
            font-style: italic;
        }

        /* 表格 */
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 1.5vw 0;
            font-size: 1.3vw;
            box-shadow: 0 2px 8px var(--shadow-color);
            border-radius: 8px;
            overflow: hidden;
        }

        th, td {
            padding: 1vw;
            text-align: left;
            border-bottom: 1px solid var(--border-color);
        }

        th {
            background: var(--button-bg);
            color: var(--header-color);
            font-weight: 600;
            text-transform: uppercase;
            font-size: 1.2vw;
            letter-spacing: 0.05em;
        }

        tr:hover {
            background: var(--button-bg);
            transition: background 0.3s ease;
        }

        /* 代码块 */
        pre, code {
            font-family: 'SF Mono', 'Monaco', 'Inconsolata', 'Fira Code', monospace;
            background: var(--code-bg);
            color: var(--code-color);
            border-radius: 6px;
        }

        code {
            padding: 0.2vw 0.5vw;
            font-size: 1.3vw;
        }

        pre {
            padding: 1.5vw;
            margin: 1.5vw 0;
            overflow-x: auto;
            font-size: 1.2vw;
            line-height: 1.6;
            box-shadow: inset 0 2px 4px var(--shadow-color);
        }

        /* 引用块 */
        blockquote {
            border-left: 4px solid var(--accent-color);
            padding: 1vw 1.5vw;
            margin: 1.5vw 0;
            background: var(--note-bg);
            border-radius: 0 8px 8px 0;
            font-style: italic;
        }

        blockquote p {
            margin: 0;
        }

        /* 提示框 */
        .note, .tip, .warning, .important {
            padding: 1.2vw;
            margin: 1.5vw 0;
            border-radius: 8px;
            border-left: 4px solid;
        }

        .note {
            background: var(--note-bg);
            border-color: var(--note-border);
        }

        .tip {
            background: #f0fdf4;
            border-color: var(--success-color);
        }

        .warning {
            background: #fffbeb;
            border-color: var(--warning-color);
        }

        .important {
            background: #fef2f2;
            border-color: var(--error-color);
        }

        /* 高亮文本 */
        .highlight {
            background: var(--highlight-bg);
            padding: 0.2vw 0.4vw;
            border-radius: 4px;
        }

        /* 卡片布局 */
        .cards {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 1.5vw;
            margin: 1.5vw 0;
        }

        .card {
            background: var(--secondary-bg);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 1.5vw;
            box-shadow: 0 2px 8px var(--shadow-color);
            transition: transform 0.3s ease, box-shadow 0.3s ease;
        }

        .card:hover {
            transform: translateY(-2px);
            box-shadow: 0 4px 16px var(--shadow-color);
        }

        .card h4 {
            margin-bottom: 0.8vw;
            color: var(--header-color);
        }

        /* 时间线 */
        .timeline {
            position: relative;
            padding-left: 3vw;
            margin: 2vw 0;
        }

        .timeline::before {
            content: '';
            position: absolute;
            left: 1vw;
            top: 0;
            bottom: 0;
            width: 2px;
            background: var(--accent-color);
        }

        .timeline-item {
            position: relative;
            margin-bottom: 2vw;
        }

        .timeline-item::before {
            content: '';
            position: absolute;
            left: -2.5vw;
            top: 0.5vw;
            width: 1vw;
            height: 1vw;
            border-radius: 50%;
            background: var(--accent-color);
            border: 3px solid var(--secondary-bg);
        }

        .timeline-date {
            font-size: 1.2vw;
            color: var(--accent-color);
            font-weight: 600;
            margin-bottom: 0.5vw;
        }

        /* 进度条 */
        .progress-bar {
            width: 100%;
            height: 2vw;
            background: var(--button-bg);
            border-radius: 1vw;
            overflow: hidden;
            margin: 1vw 0;
        }

        .progress-fill {
            height: 100%;
            background: linear-gradient(90deg, var(--accent-color), var(--header-color));
            border-radius: 1vw;
            transition: width 1s ease;
            display: flex;
            align-items: center;
            justify-content: flex-end;
            padding-right: 1vw;
            color: white;
            font-size: 1.2vw;
            font-weight: 600;
        }

        /* ==================== 控制面板 ==================== */
        #controls {
            position: fixed;
            bottom: 2vh;
            left: 50%;
            transform: translateX(-50%);
            background: var(--secondary-bg);
            backdrop-filter: blur(10px);
            padding: 1vh 2vh;
            border-radius: 50px;
            display: flex;
            align-items: center;
            gap: 1.5vh;
            z-index: 1000;
            box-shadow: 0 4px 20px var(--shadow-color);
            border: 1px solid var(--border-color);
        }

        #controls button {
            background: var(--button-bg);
            color: var(--text-color);
            border: none;
            padding: 0.8vh 1.2vh;
            border-radius: 8px;
            cursor: pointer;
            font-size: 1.8vh;
            transition: all 0.3s ease;
            display: flex;
            align-items: center;
            justify-content: center;
            min-width: 3.5vh;
            height: 3.5vh;
        }

        #controls button:hover {
            background: var(--button-hover-bg);
            transform: scale(1.1);
        }

        #controls button:active {
            transform: scale(0.95);
        }

        #slide-counter {
            font-size: 1.6vh;
            font-weight: 500;
            min-width: 8vh;
            text-align: center;
            color: var(--text-color);
            opacity: 0.8;
        }

        #progress-indicator {
            width: 20vh;
            height: 0.5vh;
            background: var(--button-bg);
            border-radius: 0.25vh;
            overflow: hidden;
        }

        #progress-bar {
            height: 100%;
            background: linear-gradient(90deg, var(--accent-color), var(--header-color));
            transition: width 0.5s ease;
        }

        /* 主题切换按钮 */
        #theme-btn {
            position: relative;
        }

        .theme-menu {
            position: absolute;
            bottom: 100%;
            left: 50%;
            transform: translateX(-50%);
            background: var(--secondary-bg);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 0.5vh;
            margin-bottom: 0.5vh;
            display: none;
            box-shadow: 0 4px 12px var(--shadow-color);
        }

        .theme-menu.show {
            display: block;
        }

        .theme-option {
            padding: 0.5vh 1vh;
            cursor: pointer;
            border-radius: 4px;
            font-size: 1.4vh;
            white-space: nowrap;
        }

        .theme-option:hover {
            background: var(--button-bg);
        }

        /* ==================== 响应式设计 ==================== */
        @media (max-width: 768px) {
            h1 { font-size: 6vw; }
            h2 { font-size: 5vw; }
            h3 { font-size: 4vw; }
            p, li { font-size: 3.5vw; }
            
            .columns, .three-columns {
                flex-direction: column;
            }
            
            .title-slide h1 { font-size: 7vw; }
            .title-slide .subtitle { font-size: 4vw; }
            
            table { font-size: 2.8vw; }
            
            #controls {
                padding: 0.5vh 1vh;
                gap: 1vh;
            }
            
            #controls button {
                padding: 0.5vh 0.8vh;
                font-size: 1.5vh;
            }
            
            #progress-indicator {
                width: 15vh;
            }
        }

        @media (min-width: 1920px) {
            h1 { font-size: 48px; }
            h2 { font-size: 36px; }
            h3 { font-size: 28px; }
            p, li { font-size: 20px; }
            
            .title-slide h1 { font-size: 64px; }
            .title-slide .subtitle { font-size: 28px; }
        }

        /* ==================== 打印样式 ==================== */
        @media print {
            body {
                background: white;
            }
            
            #presentation-container {
                box-shadow: none;
                border: none;
            }
            
            .slide {
                page-break-after: always;
                opacity: 1 !important;
                display: block !important;
                position: relative !important;
            }
            
            #controls {
                display: none;
            }
        }

        /* ==================== 动画效果 ==================== */
        @keyframes fadeIn {
            from { opacity: 0; transform: translateY(20px); }
            to { opacity: 1; transform: translateY(0); }
        }

        @keyframes slideInLeft {
            from { opacity: 0; transform: translateX(-30px); }
            to { opacity: 1; transform: translateX(0); }
        }

        @keyframes slideInRight {
            from { opacity: 0; transform: translateX(30px); }
            to { opacity: 1; transform: translateX(0); }
        }

        .fade-in {
            animation: fadeIn 0.8s ease forwards;
        }

        .slide-in-left {
            animation: slideInLeft 0.8s ease forwards;
        }

        .slide-in-right {
            animation: slideInRight 0.8s ease forwards;
        }

        /* 延迟动画 */
        .delay-1 { animation-delay: 0.2s; }
        .delay-2 { animation-delay: 0.4s; }
        .delay-3 { animation-delay: 0.6s; }
    </style>
</head>

<body>
    <div id="presentation-container">
        
        <!-- 幻灯片 1: 标题页示例 -->
        <div class="slide active">
            <div class="slide-content title-slide">
                <h1 class="fade-in">您的演示标题</h1>
                <p class="subtitle fade-in delay-1">副标题或简短描述</p>
                <p class="author fade-in delay-2"><strong>演讲者姓名</strong></p>
                <p class="institution fade-in delay-2">所属机构/部门</p>
                <p class="date fade-in delay-3">2024年12月</p>
            </div>
        </div>

        <!-- 幻灯片 2: 目录/大纲 -->
        <div class="slide">
            <div class="slide-content">
                <h1>目录</h1>
                <div class="timeline">
                    <div class="timeline-item fade-in">
                        <div class="timeline-date">第一部分</div>
                        <h3>研究背景与动机</h3>
                        <p>介绍研究领域的现状和重要性</p>
                    </div>
                    <div class="timeline-item fade-in delay-1">
                        <div class="timeline-date">第二部分</div>
                        <h3>理论基础与方法</h3>
                        <p>核心概念、理论框架和研究方法</p>
                    </div>
                    <div class="timeline-item fade-in delay-2">
                        <div class="timeline-date">第三部分</div>
                        <h3>实验结果与分析</h3>
                        <p>数据展示、结果解释和讨论</p>
                    </div>
                    <div class="timeline-item fade-in delay-3">
                        <div class="timeline-date">第四部分</div>
                        <h3>结论与展望</h3>
                        <p>主要发现、贡献和未来工作</p>
                    </div>
                </div>
            </div>
        </div>

        <!-- 幻灯片 3: 两列布局示例 -->
        <div class="slide">
            <div class="slide-content">
                <h1>两列布局示例</h1>
                <div class="columns">
                    <div class="column slide-in-left">
                        <h3>左侧内容</h3>
                        <p>这是左侧列的内容。您可以在这里放置文字、列表或其他元素。</p>
                        <ul>
                            <li>要点一：简明扼要的描述</li>
                            <li>要点二：重要信息的突出</li>
                            <li>要点三：逻辑清晰的表达</li>
                        </ul>
                        <div class="note">
                            <strong>提示：</strong>使用两列布局可以更好地组织内容，使信息呈现更加清晰。
                        </div>
                    </div>
                    <div class="column slide-in-right">
                        <h3>右侧内容</h3>
                        <figure>
                            <img src="https://via.placeholder.com/400x300/3b82f6/ffffff?text=示例图片" alt="示例图片">
                            <figcaption>图1：这里是图片说明文字</figcaption>
                        </figure>
                        <p>右侧可以放置图片、图表或补充说明。</p>
                    </div>
                </div>
            </div>
        </div>

        <!-- 幻灯片 4: 三列布局示例 -->
        <div class="slide">
            <div class="slide-content">
                <h1>三列布局与卡片</h1>
                <p>适用于比较、分类或多个并列概念的展示</p>
                <div class="cards">
                    <div class="card fade-in">
                        <h4>🎯 概念一</h4>
                        <p>第一个重要概念的简要说明。包含核心思想和关键特征。</p>
                        <div class="progress-bar">
                            <div class="progress-fill" style="width: 75%;">75%</div>
                        </div>
                    </div>
                    <div class="card fade-in delay-1">
                        <h4>🚀 概念二</h4>
                        <p>第二个重要概念的描述。展示与第一个概念的区别和联系。</p>
                        <div class="progress-bar">
                            <div class="progress-fill" style="width: 90%;">90%</div>
                        </div>
                    </div>
                    <div class="card fade-in delay-2">
                        <h4>💡 概念三</h4>
                        <p>第三个概念的解释。可以包含更多细节或应用场景。</p>
                        <div class="progress-bar">
                            <div class="progress-fill" style="width: 60%;">60%</div>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- 幻灯片 5: 数据表格示例 -->
        <div class="slide">
            <div class="slide-content">
                <h1>数据展示</h1>
                <h2>实验结果对比</h2>
                <table>
                    <thead>
                        <tr>
                            <th>方法</th>
                            <th>准确率</th>
                            <th>召回率</th>
                            <th>F1分数</th>
                            <th>时间(ms)</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td><strong>我们的方法</strong></td>
                            <td style="color: var(--success-color);">95.2%</td>
                            <td style="color: var(--success-color);">93.8%</td>
                            <td style="color: var(--success-color);">94.5%</td>
                            <td>23</td>
                        </tr>
                        <tr>
                            <td>基准方法A</td>
                            <td>89.1%</td>
                            <td>87.3%</td>
                            <td>88.2%</td>
                            <td>45</td>
                        </tr>
                        <tr>
                            <td>基准方法B</td>
                            <td>91.5%</td>
                            <td>90.2%</td>
                            <td>90.8%</td>
                            <td>31</td>
                        </tr>
                        <tr>
                            <td>传统方法</td>
                            <td>82.3%</td>
                            <td>80.1%</td>
                            <td>81.2%</td>
                            <td>67</td>
                        </tr>
                    </tbody>
                </table>
                <div class="tip">
                    <strong>结果分析：</strong>我们的方法在所有指标上都取得了最佳性能，同时保持了较快的处理速度。
                </div>
            </div>
        </div>

        <!-- 幻灯片 6: 代码示例 -->
        <div class="slide">
            <div class="slide-content">
                <h1>代码示例</h1>
                <p>展示关键算法或实现细节</p>
                <pre><code># Python 示例代码
import numpy as np
import matplotlib.pyplot as plt

class NeuralNetwork:
    def __init__(self, layers):
        """初始化神经网络"""
        self.layers = layers
        self.weights = []
        self.biases = []
        
        # 初始化权重和偏置
        for i in range(len(layers) - 1):
            w = np.random.randn(layers[i], layers[i+1]) * 0.01
            b = np.zeros((1, layers[i+1]))
            self.weights.append(w)
            self.biases.append(b)
    
    def forward(self, X):
        """前向传播"""
        self.activations = [X]
        for w, b in zip(self.weights, self.biases):
            z = np.dot(self.activations[-1], w) + b
            a = self.sigmoid(z)
            self.activations.append(a)
        return self.activations[-1]</code></pre>
                <div class="note">
                    该代码展示了一个简单的神经网络实现，包含初始化和前向传播功能。
                </div>
            </div>
        </div>

        <!-- 幻灯片 7: 数学公式示例 -->
        <div class="slide">
            <div class="slide-content">
                <h1>数学公式展示</h1>
                <h2>重要定理与推导</h2>
                <div class="columns">
                    <div class="column">
                        <h3>基础公式</h3>
                        <p>欧拉公式：</p>
                        <p style="text-align: center; font-size: 1.8vw;">
                            $$e^{i\theta} = \cos\theta + i\sin\theta$$
                        </p>
                        <p>高斯分布：</p>
                        <p style="text-align: center; font-size: 1.8vw;">
                            $$f(x) = \frac{1}{\sigma\sqrt{2\pi}} e^{-\frac{(x-\mu)^2}{2\sigma^2}}$$
                        </p>
                    </div>
                    <div class="column">
                        <h3>推导过程</h3>
                        <p>损失函数：</p>
                        <p style="text-align: center; font-size: 1.8vw;">
                            $$L = -\frac{1}{N}\sum_{i=1}^{N}[y_i\log\hat{y}_i + (1-y_i)\log(1-\hat{y}_i)]$$
                        </p>
                        <p>梯度下降更新：</p>
                        <p style="text-align: center; font-size: 1.8vw;">
                            $$\theta_{t+1} = \theta_t - \alpha \nabla_\theta L(\theta_t)$$
                        </p>
                    </div>
                </div>
            </div>
        </div>

        <!-- 幻灯片 8: 引用和参考 -->
        <div class="slide">
            <div class="slide-content">
                <h1>引用与说明</h1>
                <h2>重要观点</h2>
                <blockquote>
                    <p>"科学的真正目的不是积累事实，而是发现规律。"</p>
                    <p style="text-align: right; font-size: 1.2vw;">—— 威廉·劳伦斯·布拉格</p>
                </blockquote>
                
                <h3>不同类型的提示框</h3>
                <div class="note">
                    <strong>📝 注意：</strong>这是一个普通的提示信息，用于补充说明。
                </div>
                <div class="tip">
                    <strong>💡 提示：</strong>这是一个有用的技巧或建议。
                </div>
                <div class="warning">
                    <strong>⚠️ 警告：</strong>这是需要特别注意的事项。
                </div>
                <div class="important">
                    <strong>❗ 重要：</strong>这是关键信息，请务必关注。
                </div>
            </div>
        </div>

        <!-- 幻灯片 9: 图文混排 -->
        <div class="slide">
            <div class="slide-content">
                <h1>研究方法</h1>
                <div class="columns">
                    <div class="column column-60">
                        <h3>实验设计</h3>
                        <ol>
                            <li><strong>数据收集阶段</strong>
                                <ul>
                                    <li>收集原始数据</li>
                                    <li>数据清洗和预处理</li>
                                    <li>特征工程</li>
                                </ul>
                            </li>
                            <li><strong>模型训练阶段</strong>
                                <ul>
                                    <li>选择合适的算法</li>
                                    <li>参数调优</li>
                                    <li>交叉验证</li>
                                </ul>
                            </li>
                            <li><strong>评估阶段</strong>
                                <ul>
                                    <li>性能评估</li>
                                    <li>结果分析</li>
                                </ul>
                            </li>
                        </ol>
                    </div>
                    <div class="column column-40">
                        <figure>
                            <img src="https://via.placeholder.com/300x400/10b981/ffffff?text=流程图" alt="实验流程">
                            <figcaption>图2：实验流程示意图</figcaption>
                        </figure>
                    </div>
                </div>
            </div>
        </div>

        <!-- 幻灯片 10: 总结页 -->
        <div class="slide">
            <div class="slide-content">
                <h1>总结与展望</h1>
                <div class="columns">
                    <div class="column">
                        <h3>🎯 主要贡献</h3>
                        <ul>
                            <li>提出了创新的解决方案</li>
                            <li>实现了性能的显著提升</li>
                            <li>验证了理论的正确性</li>
                            <li>拓展了应用领域</li>
                        </ul>
                    </div>
                    <div class="column">
                        <h3>🚀 未来工作</h3>
                        <ul>
                            <li>扩展到更多应用场景</li>
                            <li>优化算法效率</li>
                            <li>深入理论研究</li>
                            <li>开发实用工具</li>
                        </ul>
                    </div>
                </div>
                <div class="important" style="margin-top: 3vw;">
                    <strong>核心观点：</strong>本研究为该领域提供了新的视角和方法，具有重要的理论意义和实践价值。
                </div>
            </div>
        </div>

        <!-- 幻灯片 11: 致谢页 -->
        <div class="slide">
            <div class="slide-content title-slide">
                <h1>谢谢！</h1>
                <p class="subtitle">感谢您的聆听</p>
                <br>
                <h2 style="color: var(--text-color);">问答环节</h2>
                <p style="font-size: 1.5vw; margin-top: 2vw;">
                    📧 your.email@example.com<br>
                    🌐 your-website.com<br>
                    📱 微信/GitHub: YourID
                </p>
            </div>
        </div>

        <!-- 幻灯片 12: 使用说明 -->
        <div class="slide">
            <div class="slide-content">
                <h1>模板使用说明</h1>
                <h2>快捷键操作</h2>
                <div class="columns">
                    <div class="column">
                        <h3>导航控制</h3>
                        <table style="font-size: 1.2vw;">
                            <tr><td><code>→</code> / <code>Space</code></td><td>下一页</td></tr>
                            <tr><td><code>←</code></td><td>上一页</td></tr>
                            <tr><td><code>Home</code></td><td>第一页</td></tr>
                            <tr><td><code>End</code></td><td>最后一页</td></tr>
                            <tr><td><code>G</code></td><td>跳转到指定页</td></tr>
                            <tr><td><code>F</code></td><td>全屏切换</td></tr>
                        </table>
                    </div>
                    <div class="column">
                        <h3>功能特性</h3>
                        <ul>
                            <li>🎨 <strong>主题切换</strong>：点击调色板按钮切换主题</li>
                            <li>📱 <strong>响应式设计</strong>：自适应不同屏幕尺寸</li>
                            <li>🔢 <strong>数学公式</strong>：支持LaTeX数学公式渲染</li>
                            <li>📊 <strong>丰富组件</strong>：表格、代码、时间线等</li>
                            <li>🎯 <strong>动画效果</strong>：淡入、滑入等过渡动画</li>
                            <li>🖨️ <strong>打印友好</strong>：优化的打印样式</li>
                        </ul>
                    </div>
                </div>
                <div class="tip">
                    <strong>提示：</strong>您可以根据需要修改CSS变量来自定义配色方案，或添加新的幻灯片布局样式。
                </div>
            </div>
        </div>

    </div>

    <!-- 控制面板 -->
    <div id="controls">
        <button id="prev-btn" title="上一页 (←)">◀</button>
        <button id="first-btn" title="第一页 (Home)">⏮</button>
        <span id="slide-counter">1 / 12</span>
        <button id="last-btn" title="最后一页 (End)">⏭</button>
        <button id="next-btn" title="下一页 (→)">▶</button>
        <div id="progress-indicator">
            <div id="progress-bar" style="width: 8.33%;"></div>
        </div>
        <button id="goto-btn" title="跳转 (G)">🔢</button>
        <button id="theme-btn" title="切换主题">🎨
            <div class="theme-menu" id="theme-menu">
                <div class="theme-option" data-theme="default">专业蓝</div>
                <div class="theme-option" data-theme="dark">深色</div>
                <div class="theme-option" data-theme="green">清新绿</div>
            </div>
        </button>
        <button id="fullscreen-btn" title="全屏 (F)">⛶</button>
    </div>

    <!-- JavaScript 控制逻辑 -->
    <script>
        document.addEventListener('DOMContentLoaded', function() {
            // 获取元素
            const slides = document.querySelectorAll('.slide');
            const prevBtn = document.getElementById('prev-btn');
            const nextBtn = document.getElementById('next-btn');
            const firstBtn = document.getElementById('first-btn');
            const lastBtn = document.getElementById('last-btn');
            const gotoBtn = document.getElementById('goto-btn');
            const fullscreenBtn = document.getElementById('fullscreen-btn');
            const themeBtn = document.getElementById('theme-btn');
            const themeMenu = document.getElementById('theme-menu');
            const slideCounter = document.getElementById('slide-counter');
            const progressBar = document.getElementById('progress-bar');
            const presentationContainer = document.getElementById('presentation-container');
            
            let currentSlide = 0;
            const totalSlides = slides.length;
            
            // 更新幻灯片显示
            function showSlide(index) {
                if (index < 0) index = 0;
                if (index >= totalSlides) index = totalSlides - 1;
                
                slides.forEach((slide, i) => {
                    slide.classList.remove('active');
                    if (i === index) {
                        slide.classList.add('active');
                        // 触发动画
                        const animatedElements = slide.querySelectorAll('.fade-in, .slide-in-left, .slide-in-right');
                        animatedElements.forEach(el => {
                            el.style.animation = 'none';
                            el.offsetHeight; // 触发重排
                            el.style.animation = null;
                        });
                    }
                });
                
                currentSlide = index;
                updateControls();
                
                // 渲染数学公式
                if (window.MathJax && window.MathJax.typeset) {
                    MathJax.typesetPromise([slides[index]]).catch(err => {
                        console.error('MathJax error:', err);
                    });
                }
            }
            
            // 更新控制面板
            function updateControls() {
                slideCounter.textContent = `${currentSlide + 1} / ${totalSlides}`;
                progressBar.style.width = `${((currentSlide + 1) / totalSlides) * 100}%`;
                
                // 更新按钮状态
                prevBtn.disabled = currentSlide === 0;
                firstBtn.disabled = currentSlide === 0;
                nextBtn.disabled = currentSlide === totalSlides - 1;
                lastBtn.disabled = currentSlide === totalSlides - 1;
            }
            
            // 导航功能
            function nextSlide() {
                if (currentSlide < totalSlides - 1) {
                    showSlide(currentSlide + 1);
                }
            }
            
            function prevSlide() {
                if (currentSlide > 0) {
                    showSlide(currentSlide - 1);
                }
            }
            
            function firstSlide() {
                showSlide(0);
            }
            
            function lastSlide() {
                showSlide(totalSlides - 1);
            }
            
            function gotoSlide() {
                const slideNum = prompt(`跳转到第几页？(1-${totalSlides})`);
                if (slideNum) {
                    const num = parseInt(slideNum) - 1;
                    if (!isNaN(num) && num >= 0 && num < totalSlides) {
                        showSlide(num);
                    }
                }
            }
            
            // 全屏功能
            function toggleFullscreen() {
                if (!document.fullscreenElement) {
                    presentationContainer.requestFullscreen().catch(err => {
                        console.error('Fullscreen error:', err);
                    });
                } else {
                    document.exitFullscreen();
                }
            }
            
            // 主题切换
            function setTheme(theme) {
                document.documentElement.setAttribute('data-theme', theme);
                localStorage.setItem('presentation-theme', theme);
                themeMenu.classList.remove('show');
            }
            
            // 事件监听
            prevBtn.addEventListener('click', prevSlide);
            nextBtn.addEventListener('click', nextSlide);
            firstBtn.addEventListener('click', firstSlide);
            lastBtn.addEventListener('click', lastSlide);
            gotoBtn.addEventListener('click', gotoSlide);
            fullscreenBtn.addEventListener('click', toggleFullscreen);
            
            // 主题菜单
            themeBtn.addEventListener('click', (e) => {
                e.stopPropagation();
                themeMenu.classList.toggle('show');
            });
            
            document.addEventListener('click', () => {
                themeMenu.classList.remove('show');
            });
            
            document.querySelectorAll('.theme-option').forEach(option => {
                option.addEventListener('click', (e) => {
                    e.stopPropagation();
                    setTheme(option.dataset.theme);
                });
            });
            
            // 键盘事件
            document.addEventListener('keydown', (e) => {
                // 检查是否在可滚动内容中
                if (e.target.closest('.slide-content')) {
                    const content = e.target.closest('.slide-content');
                    if (content.scrollHeight > content.clientHeight) {
                        // 允许在可滚动内容中使用方向键
                        if (e.key === 'ArrowUp' || e.key === 'ArrowDown') {
                            return;
                        }
                    }
                }
                
                switch(e.key) {
                    case 'ArrowRight':
                    case ' ':
                    case 'PageDown':
                        e.preventDefault();
                        nextSlide();
                        break;
                    case 'ArrowLeft':
                    case 'PageUp':
                        e.preventDefault();
                        prevSlide();
                        break;
                    case 'Home':
                        e.preventDefault();
                        firstSlide();
                        break;
                    case 'End':
                        e.preventDefault();
                        lastSlide();
                        break;
                    case 'f':
                    case 'F':
                        e.preventDefault();
                        toggleFullscreen();
                        break;
                    case 'g':
                    case 'G':
                        e.preventDefault();
                        gotoSlide();
                        break;
                    case 'Escape':
                        if (document.fullscreenElement) {
                            document.exitFullscreen();
                        }
                        break;
                }
            });
            
            // 触摸手势支持
            let touchStartX = 0;
            let touchEndX = 0;
            
            presentationContainer.addEventListener('touchstart', (e) => {
                touchStartX = e.changedTouches[0].screenX;
            });
            
            presentationContainer.addEventListener('touchend', (e) => {
                touchEndX = e.changedTouches[0].screenX;
                handleSwipe();
            });
            
            function handleSwipe() {
                const swipeThreshold = 50;
                const diff = touchStartX - touchEndX;
                
                if (Math.abs(diff) > swipeThreshold) {
                    if (diff > 0) {
                        nextSlide(); // 左滑
                    } else {
                        prevSlide(); // 右滑
                    }
                }
            }
            
            // 鼠标滚轮支持（可选）
            let scrollTimeout;
            presentationContainer.addEventListener('wheel', (e) => {
                // 检查是否在可滚动内容中
                const content = e.target.closest('.slide-content');
                if (content && content.scrollHeight > content.clientHeight) {
                    // 允许在内容中滚动
                    return;
                }
                
                e.preventDefault();
                clearTimeout(scrollTimeout);
                scrollTimeout = setTimeout(() => {
                    if (e.deltaY > 0) {
                        nextSlide();
                    } else {
                        prevSlide();
                    }
                }, 50);
            }, { passive: false });
            
            // 加载保存的主题
            const savedTheme = localStorage.getItem('presentation-theme');
            if (savedTheme) {
                setTheme(savedTheme);
            }
            
            // 初始化
            showSlide(0);
            
            // URL参数支持（可以通过URL指定起始页）
            const urlParams = new URLSearchParams(window.location.search);
            const startSlide = urlParams.get('slide');
            if (startSlide) {
                const slideNum = parseInt(startSlide) - 1;
                if (!isNaN(slideNum) && slideNum >= 0 && slideNum < totalSlides) {
                    showSlide(slideNum);
                }
            }
            
            // 自动播放功能（可选）
            let autoplayInterval;
            let isAutoplay = false;
            
            function toggleAutoplay() {
                if (isAutoplay) {
                    clearInterval(autoplayInterval);
                    isAutoplay = false;
                } else {
                    autoplayInterval = setInterval(() => {
                        if (currentSlide < totalSlides - 1) {
                            nextSlide();
                        } else {
                            clearInterval(autoplayInterval);
                            isAutoplay = false;
                        }
                    }, 5000); // 5秒切换
                    isAutoplay = true;
                }
            }
            
            // 可以添加自动播放按钮
            // document.addEventListener('keydown', (e) => {
            //     if (e.key === 'a' || e.key === 'A') {
            //         toggleAutoplay();
            //     }
            // });
        });
    </script>
</body>

</html>