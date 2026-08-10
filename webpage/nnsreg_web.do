*! nnsreg web page and interactive companion
clear all
set more off
set scheme s1color

local root "/Users/ericbooth/Library/CloudStorage/GoogleDrive-eric.booth@texas2036.org/My Drive/NNS-stata-public"
cd "`root'"
adopath + "code"

cap mkdir "webpage/assets"

// Build a statashiny companion. This is an explanatory browser-side
// visualization, not a separate Stata estimator.
statashiny, title("NNS Fit Explorer") ///
    subtitle("How partition detail changes a line-segment fit") replace
statashiny order, input(num) label("Partition order") val("3") min(1) max(6) step(1)
statashiny threshold, input(num) label("Threshold") val("35") min(15) max(60) step(1)
statashiny contrast, input(num) label("Post-threshold slope") val("0.55") min(0.1) max(1.2) step(0.05)
statashiny roughness, input(num) label("Noise amplitude") val("2") min(0) max(6) step(0.5)
statashiny, output(raw) label("<style>#nnsCanvas{width:100%;height:360px;display:block}.fit-card{max-width:100%;overflow:hidden}.fit-note{max-width:68ch}</style><div class='card statashiny-component'><div class='card-body'><h6 class='card-title'>How to read this explorer</h6><p class='small text-muted mb-1'>A teaching visualization from the <em>nnsreg</em> Stata Journal article: how partition <strong>order</strong> trades resolution against noise. It does not run <code>nnsreg</code>.</p><p class='small text-muted mb-1'><strong>Sliders:</strong> Threshold, Post-threshold slope, and Noise amplitude shape the true (navy) curve; Partition order sets how many line segments the (orange) fit uses.</p><p class='small text-muted mb-1'><strong>Try:</strong> raise the order to hug the curve; then raise noise and order together to watch the fit chase noise (overfitting).</p><p class='small text-muted mb-0'><strong>Caveat:</strong> here the segments are evenly spaced and drawn on the curve; the real <code>nnsreg</code> places its knots adaptively from the data and fits the observations, so its knots concentrate where the data bend rather than at even intervals.</p></div></div><div class='card statashiny-component fit-card'><div class='card-body'><h5 class='card-title'>Interactive fit</h5><canvas id='nnsCanvas'></canvas><p id='fitNote' class='text-muted small mt-3 fit-note'></p></div></div>")
statashiny, calc("const canvas=document.getElementById('nnsCanvas'); if(!canvas) return; const ctx=canvas.getContext('2d'); const parent=canvas.parentElement; const W=Math.max(520,parent.clientWidth-8); const H=360; canvas.width=W; canvas.height=H; const order=Math.max(1,Math.min(6,parseInt(document.getElementById('order').value||3))); const threshold=+document.getElementById('threshold').value; const contrast=+document.getElementById('contrast').value; const rough=+document.getElementById('roughness').value; function f(x){const base=8+0.02*x; const jump=x>threshold?contrast*(x-threshold)+0.006*Math.pow(x-threshold,2):0; return base+jump+rough*Math.sin(x/4);} const margin={l:52,r:24,t:22,b:42}; function sx(x){return margin.l+(x/80)*(W-margin.l-margin.r);} function sy(y){return H-margin.b-(y/55)*(H-margin.t-margin.b);} ctx.clearRect(0,0,W,H); ctx.fillStyle='#ffffff'; ctx.fillRect(0,0,W,H); ctx.strokeStyle='#d8dee8'; ctx.lineWidth=1; ctx.font='12px sans-serif'; ctx.fillStyle='#364152'; for(let gx=0;gx<=80;gx+=20){ctx.beginPath();ctx.moveTo(sx(gx),margin.t);ctx.lineTo(sx(gx),H-margin.b);ctx.stroke();ctx.fillText(gx,sx(gx)-6,H-18);} for(let gy=0;gy<=50;gy+=10){ctx.beginPath();ctx.moveTo(margin.l,sy(gy));ctx.lineTo(W-margin.r,sy(gy));ctx.stroke();ctx.fillText(gy,14,sy(gy)+4);} ctx.strokeStyle='#6c7a8d'; ctx.lineWidth=1.5; ctx.beginPath(); ctx.moveTo(margin.l,margin.t); ctx.lineTo(margin.l,H-margin.b); ctx.lineTo(W-margin.r,H-margin.b); ctx.stroke(); ctx.fillStyle='#6c7a8d'; for(let i=0;i<90;i++){const x=(i*37%80)+((i%5)-2)*0.22; const y=f(x)+((i*19)%9-4)*0.55; ctx.beginPath(); ctx.arc(sx(x),sy(y),3,0,Math.PI*2); ctx.strokeStyle='#6c7a8d'; ctx.stroke();} ctx.strokeStyle='#1b2d55'; ctx.lineWidth=2; ctx.beginPath(); for(let x=0;x<=80;x+=1){const y=f(x); if(x===0) ctx.moveTo(sx(x),sy(y)); else ctx.lineTo(sx(x),sy(y));} ctx.stroke(); const pieces=Math.pow(2,order-1); ctx.strokeStyle='#d44500'; ctx.lineWidth=4; ctx.beginPath(); for(let j=0;j<=pieces;j++){const x=80*j/pieces; const y=f(x); if(j===0) ctx.moveTo(sx(x),sy(y)); else ctx.lineTo(sx(x),sy(y));} ctx.stroke(); ctx.fillStyle='#1b2d55'; ctx.fillText('Smooth reference curve', W-180, 24); ctx.fillStyle='#d44500'; ctx.fillText('Line-segment fit', W-180, 44); document.getElementById('fitNote').innerText='Order '+order+' uses '+pieces+' line segments in this teaching display. Higher order follows local shape more closely, but it can also track noise.';")
statashiny, build export("webpage/nnsreg_interactive.html")

cd "webpage"
wdinit nnsreg_web, replace

webdoc put <h1>nnsreg: Reportable segment slopes for nonlinear relationships in Stata</h1>
webdoc put <p><strong>nnsreg</strong> fits a univariate NNS-style piecewise linear curve and posts the result as a Stata estimator. Use it to inspect bends, flat regions, slope changes, and whether those features survive sensitivity checks. The examples here match the help file, README, and Stata Journal draft.</p>
webdoc put <p>This page is a companion to the Stata Journal article <em>nnsreg: Reportable segment slopes for nonlinear relationships in Stata</em>. The article makes a simple point: a flexible smoother can draw a nonlinear curve, but applied work often needs a few reportable statements about where a relationship is flat, steep, or changing, and <code>nnsreg</code> returns exactly that as ordinary Stata estimation results. The sections below follow the same path as the article: the questions the command answers, its syntax, the simulated designs that let us score it against a known curve, three worked examples, sensitivity diagnostics, and an interactive explorer for intuition. In each worked example, read the estimates table for the named first-segment slope and the slope-change coefficients, and read the figure for how well each method tracks the true curve. A kernel smoother often recovers a smooth curve at least as accurately; the point of <code>nnsreg</code> is the reportable segment summary it hands back, not a better fit.</p>

webdoc put <h2>Questions It Helps Answer</h2>
webdoc put <p>Use <code>nnsreg</code> when you need to ask where a relationship bends, which part of the support is steep or flat, whether slope changes are stable across options, and whether a segmented fit adds value relative to <code>regress</code> or <code>npregress kernel</code>. Kernel smoothing remains a strong choice for smooth conditional-mean recovery; NNS is most useful when segmentation, clustered support, sharp local features, or stored slope-change output matter.</p>

webdoc put <h2>Syntax And Options</h2>
webdoc put <pre><code>nnsreg depvar indepvar [if] [in] [, order(integer) noise(mean|median) method(connect|ols) partition(xy|xonly) minobs(integer) level(cilevel) vce(robust|cluster clustvar) slopes noplot knotlines graph_opts(string) generate(newvar)]</code></pre>
webdoc put <ul><li><code>order(#)</code> controls partition depth. Larger values allow more line segments.</li><li><code>noise(mean)</code> uses means for centers and regression points.</li><li><code>noise(median)</code> is a robustness check for skewness or outliers.</li><li><code>method(connect)</code> connects NNS regression points directly.</li><li><code>method(ols)</code> estimates an OLS spline at NNS-derived knots.</li><li><code>partition(xy)</code> is the default NNS-style quadrant partition.</li><li><code>partition(xonly)</code> is a vertical-band sensitivity check.</li><li><code>minobs(#)</code> sets the smallest group the partition keeps splitting (default 8).</li><li><code>slopes</code> reports each segment's slope and standard error, and stores them in <code>e(segments)</code>.</li><li><code>vce(robust|cluster clustvar)</code> passes a robust or cluster-robust covariance to the least-squares stage (with <code>method(ols)</code>).</li></ul>
webdoc put <p>The reporting payoff lives in the stored results: <code>e(knots)</code> holds the data-derived knot locations, <code>e(b)</code> names the first-segment slope and each slope change, and <code>e(segments)</code> (or the <code>slopes</code> option) reports the cumulative slope in each range. These are what let you turn a fitted curve into a few <code>lincom</code>-testable statements, the article's central contribution. One honesty note carried over from the article: for <code>method(connect)</code> the standard errors are heuristic references from a rescaled OLS spline covariance, so read them as descriptive summaries of segment stability rather than as tests.</p>

webdoc put <h2>Simulation Design</h2>
webdoc put <p>The education example simulates an S-shaped curve with diminishing returns. The water example simulates a flatter lower range and a sharper increase after an infrastructure-age threshold. The clustered peak/dip design simulates dense local regions with a peak, sharp dip, and recovery. These are simulated examples, so the figures test software behavior against known curves rather than estimating real policy effects.</p>
webdoc put <p>Why simulate? Because the true curve is known in each design, the figures can score how well a method recovers the signal, not just how closely it fits the noisy points. That is the distinction the article measures with its true-curve RMSE. Each design also isolates a different challenge: a smooth curve (education), a threshold or kink (water), and clustered support with a sharp local feature (the clustered design). Watch how the same command behaves across all three.</p>

webdoc put <h2>Worked Example: Education</h2>
button
cd "`root'"
adopath + "code"
use "examples/texas_education.dta", clear
regress staar_pass spend_pupil
npregress kernel staar_pass spend_pupil
nnsreg staar_pass spend_pupil, order(3) method(connect) partition(xy) generate(fit_connect)
estimates store edu_connect
nnsreg staar_pass spend_pupil, order(3) method(ols) partition(xy) generate(fit_ols)
estimates store edu_ols
estimates table edu_connect edu_ols, b se stats(r2 rmse N)
buttonclose
webdoc put <figure><img src="../manuscript/fig_edu_fits.png" class="img-fluid" alt="Education fit comparison, one method per panel"><figcaption>One fitted method per panel against the same simulated data (gray) and true curve (black). The kernel smoother and NNS OLS both recover the S-shape; the connected NNS fit is more jagged because it interpolates the regression points.</figcaption></figure>
webdoc put <figure><img src="../manuscript/fig_edu_coef.png" class="img-fluid" alt="Education stored slope components"><figcaption>Stored slope components with knot-located labels: the reportable first-segment slope and slope changes that <code>e(b)</code>, <code>e(segments)</code>, and <code>lincom</code> return.</figcaption></figure>
webdoc put <p><strong>Why this example:</strong> the education design is a smooth S-curve, the case where a kernel smoother should do at least as well as <code>nnsreg</code>. It shows what the command adds even when it does not win on accuracy. <strong>How to read the output:</strong> in the estimates table, the first coefficient (<code>spend_pupil</code>) is the slope of the first segment, and each <code>knot#</code> coefficient is the change in slope at an NNS-derived knot; their locations are in <code>e(knots)</code>. In the figure, the kernel smoother and NNS OLS both recover the S-shape, while the connected NNS fit is more jagged because it interpolates the regression points rather than minimizing squared error. The lesson matches the article: on a smooth curve the smoother matches or beats NNS on accuracy, and what NNS adds is the reportable slope pieces, not a better fit.</p>

webdoc put <h2>Worked Example: Water</h2>
button
use "examples/texas_water.dta", clear
regress water_loss pipe_age
npregress kernel water_loss pipe_age
nnsreg water_loss pipe_age, order(3) method(connect) partition(xy) generate(fit_connect)
estimates store water_connect
nnsreg water_loss pipe_age, order(3) method(ols) partition(xy) generate(fit_ols)
estimates store water_ols
estimates table water_connect water_ols, b se stats(r2 rmse N)
buttonclose
webdoc put <figure><img src="../manuscript/fig_water_fits.png" class="img-fluid" alt="Water fit comparison, one method per panel"><figcaption>One fitted method per panel. The true curve is flat until pipe age 35 and rises afterward; the kernel smoother and NNS OLS both recover the flat-then-rising shape.</figcaption></figure>
webdoc put <figure><img src="../manuscript/fig_water_knots.png" class="img-fluid" alt="Water NNS OLS fit with interior knots"><figcaption>The NNS OLS fit with its interior knots (from <code>e(knots)</code>, dashed): several land near the true change point at age 35.</figcaption></figure>
webdoc put <p><strong>Why this example:</strong> the water curve is flat for young pipes and rises after a change point near age 35, so it tests whether the command can locate a threshold. <strong>How to read the output:</strong> read <code>e(knots)</code> and the figure to see several NNS knots land near that change point; the first-segment slope is close to zero (flat young pipes) and the terminal slope is clearly positive (deteriorating old pipes). That flat-then-rising contrast is a single <code>lincom</code> on the stored coefficients. This is the interpretive case the article highlights: the command reports <em>where</em> the relationship turns and how steep it is on each side, not just that it curves.</p>

webdoc put <h2>Worked Example: Clustered Peak/Dip</h2>
button
use "examples/clustered_peak_dip.dta", clear
regress outcome_index cluster_score
npregress kernel outcome_index cluster_score
nnsreg outcome_index cluster_score, order(3) method(connect) partition(xy) generate(fit_connect)
estimates store cluster_connect
nnsreg outcome_index cluster_score, order(4) method(ols) partition(xy) generate(fit_ols)
estimates store cluster_ols
estimates table cluster_connect cluster_ols, b se stats(r2 rmse N)
buttonclose
webdoc put <figure><img src="../manuscript/fig_cluster_fits.png" class="img-fluid" alt="Clustered peak/dip fit comparison, one method per panel"><figcaption>One fitted method per panel on the clustered peak/dip design, where NNS OLS improves recovery of the known curve relative to a single kernel bandwidth.</figcaption></figure>
webdoc put <p><strong>Why this example:</strong> observations cluster in three regions with a sharp dip between them, the setting the NNS literature identifies as favorable to partition-based fits, because a single kernel bandwidth must compromise across regions of very different local density. <strong>How to read the output:</strong> NNS OLS recovers the known curve more closely than the kernel here. The larger point comes from the article's 100-replication Monte Carlo: in this design the segmented fit matches the best smoother's average accuracy while varying far less from sample to sample, so its advantage is reliability rather than a lower average error. This example uses <code>order(4)</code>, one step finer than the default, because the sharp local features reward a little more resolution.</p>

webdoc put <h2>Diagnostics And Sensitivity</h2>
webdoc put <p>The verification script exports <code>table_model_diagnostics.csv</code> with a common in-sample <code>R^2</code>, observed RMSE, and true-curve RMSE. In the current run, <code>npregress kernel</code> and NNS OLS are close in the smooth education and water simulations, while NNS OLS improves substantially in the clustered peak/dip design.</p>
webdoc put <p><strong>How to read the sensitivity figures:</strong> each one overlays fits across <code>order()</code> and <code>noise()</code> choices. The question to ask is whether the substantive shape survives the variations: the plateau in education, the flat-then-rising threshold in water, the peak and dip in the clustered design. When a feature appears in every variant, it is safe to report; when it disappears between <code>order(2)</code> and <code>order(3)</code>, or between mean and median centers, report it as tentative. The clustered figure also shows the failure mode: at <code>order(5)</code> the fit adds local movement that improves in-sample fit but drifts from the true curve, the overfitting the article warns against. Running this small grid before reporting a bend is the habit the article recommends.</p>
webdoc put <figure><img src="../manuscript/fig_edu_sensitivity.png" class="img-fluid" alt="Education sensitivity figure"><figcaption>Education sensitivity to <code>order()</code> and <code>noise()</code>.</figcaption></figure>
webdoc put <figure><img src="../manuscript/fig_water_sensitivity.png" class="img-fluid" alt="Water sensitivity figure"><figcaption>Water sensitivity to <code>order()</code> and <code>noise()</code>.</figcaption></figure>
webdoc put <figure><img src="../manuscript/fig_cluster_sensitivity.png" class="img-fluid" alt="Clustered peak/dip sensitivity figure"><figcaption>Clustered peak/dip sensitivity to <code>order()</code>.</figcaption></figure>

webdoc put <h2>Real-Data Illustration: Motorcycle Impact</h2>
webdoc put <p>The simulated designs use known curves so each method can be scored against the truth. To watch <code>nnsreg</code> on data it was not built around, the article turns to the motorcycle-impact benchmark: accelerometer readings against time after impact (<code>webuse motorcycle</code>). The regressor is heavily tied and clustered, and the response has a flat run, a sharp trough, and a rebound, so this is the command's favorable case on real data.</p>
webdoc put <figure><img src="../manuscript/fig_moto_fits.png" class="img-fluid" alt="Motorcycle fit comparison, one method per panel"><figcaption>One fitted method per panel against the 133 accelerometer readings. The kernel smoother and NNS OLS both recover the trough and rebound; the connected NNS fit is rougher.</figcaption></figure>
webdoc put <figure><img src="../manuscript/fig_moto_knots.png" class="img-fluid" alt="Motorcycle NNS OLS fit with interior knots"><figcaption>The NNS OLS fit with its interior knots (from <code>e(knots)</code>), which concentrate in the impact trough where the response changes fastest.</figcaption></figure>
webdoc put <p>Because the true curve is unknown here, the article scores accuracy by repeated cross-validation. On held-out error the NNS OLS spline and the kernel smoother finish close together, with <code>npregress series</code> lowest; the value of <code>nnsreg</code> is the reportable segment slopes it hands back and its ability to score held-out points directly from stored results. Reproduce it with <code>motorcycle_example.do</code>.</p>

webdoc put <h2>Interactive Fit Explorer</h2>
webdoc put <p>The explorer below is a teaching visualization for one idea from the article: the tradeoff between resolution and noise that the <code>order()</code> option controls. It does not run <code>nnsreg</code> and it is not a Stata estimation result. Instead it lets you set a known curve and watch how a coarser or finer piecewise-linear fit tracks it, so the sensitivity figures above become intuitive.</p>
webdoc put <p><strong>What the sliders do.</strong> Three of them shape the <em>true</em> curve you are trying to recover: <strong>Threshold</strong> sets where the curve kinks (as pipe age 35 does in the water example), <strong>Post-threshold slope</strong> sets how steeply it rises after the kink, and <strong>Noise amplitude</strong> adds wiggle around it. The fourth, <strong>Partition order</strong>, controls the <em>fit</em>: raising it lets the line-segment fit use more pieces.</p>
webdoc put <p><strong>How to read the canvas.</strong> The gray dots are synthetic observations, the navy line is the smooth reference curve you set with the sliders, and the orange line is the piecewise-linear fit. The note under the canvas reports how many segments the current order uses.</p>
webdoc put <p><strong>What to try.</strong> Raise <strong>Partition order</strong> from 1 toward 6 and watch the orange fit hug the navy curve more closely, since more segments follow local shape. Then raise <strong>Noise amplitude</strong> and keep the order high, and watch the orange fit start chasing the noise instead of the curve. That is overfitting, the same failure the article shows at <code>order(5)</code> in the clustered design, and the reason the default is <code>order(3)</code> paired with a sensitivity check rather than the highest order available.</p>
webdoc put <p><strong>What this display simplifies.</strong> To stay fast in a browser, the segments here are placed at even spacing along the curve, and the fit is drawn on the reference curve rather than estimated from the scattered dots. The real <code>nnsreg</code> differs in exactly the way that matters: it places its knots <em>adaptively</em> from the data through the partial-moment partition, so they concentrate where the data bend rather than at even intervals, and it fits the actual observations. Use the worked examples and figures above for how the command truly behaves; use this explorer only to feel the detail-versus-noise tradeoff.</p>
webdoc put <div style="max-width:100%;overflow:hidden"><iframe src="nnsreg_interactive.html" style="width:100%;height:900px;border:1px solid #d8dee8;border-radius:8px;"></iframe></div>

webdoc put <h2>Install</h2>
webdoc put <pre><code>net install nnsreg, from("https://raw.githubusercontent.com/ericabooth/NNS-stata-public/main/github") replace</code></pre>

webdoc close
