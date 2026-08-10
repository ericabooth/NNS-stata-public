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
statashiny, output(raw) label("<style>#nnsCanvas{width:100%;height:360px;display:block}.fit-card{max-width:100%;overflow:hidden}.fit-note{max-width:68ch}</style><div class='card statashiny-component fit-card'><div class='card-body'><h5 class='card-title'>Interactive fit</h5><canvas id='nnsCanvas'></canvas><p id='fitNote' class='text-muted small mt-3 fit-note'></p></div></div>")
statashiny, calc("const canvas=document.getElementById('nnsCanvas'); if(!canvas) return; const ctx=canvas.getContext('2d'); const parent=canvas.parentElement; const W=Math.max(520,parent.clientWidth-8); const H=360; canvas.width=W; canvas.height=H; const order=Math.max(1,Math.min(6,parseInt(document.getElementById('order').value||3))); const threshold=+document.getElementById('threshold').value; const contrast=+document.getElementById('contrast').value; const rough=+document.getElementById('roughness').value; function f(x){const base=8+0.02*x; const jump=x>threshold?contrast*(x-threshold)+0.006*Math.pow(x-threshold,2):0; return base+jump+rough*Math.sin(x/4);} const margin={l:52,r:24,t:22,b:42}; function sx(x){return margin.l+(x/80)*(W-margin.l-margin.r);} function sy(y){return H-margin.b-(y/55)*(H-margin.t-margin.b);} ctx.clearRect(0,0,W,H); ctx.fillStyle='#ffffff'; ctx.fillRect(0,0,W,H); ctx.strokeStyle='#d8dee8'; ctx.lineWidth=1; ctx.font='12px sans-serif'; ctx.fillStyle='#364152'; for(let gx=0;gx<=80;gx+=20){ctx.beginPath();ctx.moveTo(sx(gx),margin.t);ctx.lineTo(sx(gx),H-margin.b);ctx.stroke();ctx.fillText(gx,sx(gx)-6,H-18);} for(let gy=0;gy<=50;gy+=10){ctx.beginPath();ctx.moveTo(margin.l,sy(gy));ctx.lineTo(W-margin.r,sy(gy));ctx.stroke();ctx.fillText(gy,14,sy(gy)+4);} ctx.strokeStyle='#6c7a8d'; ctx.lineWidth=1.5; ctx.beginPath(); ctx.moveTo(margin.l,margin.t); ctx.lineTo(margin.l,H-margin.b); ctx.lineTo(W-margin.r,H-margin.b); ctx.stroke(); ctx.fillStyle='#6c7a8d'; for(let i=0;i<90;i++){const x=(i*37%80)+((i%5)-2)*0.22; const y=f(x)+((i*19)%9-4)*0.55; ctx.beginPath(); ctx.arc(sx(x),sy(y),3,0,Math.PI*2); ctx.strokeStyle='#6c7a8d'; ctx.stroke();} ctx.strokeStyle='#1b2d55'; ctx.lineWidth=2; ctx.beginPath(); for(let x=0;x<=80;x+=1){const y=f(x); if(x===0) ctx.moveTo(sx(x),sy(y)); else ctx.lineTo(sx(x),sy(y));} ctx.stroke(); const pieces=Math.pow(2,order-1); ctx.strokeStyle='#d44500'; ctx.lineWidth=4; ctx.beginPath(); for(let j=0;j<=pieces;j++){const x=80*j/pieces; const y=f(x); if(j===0) ctx.moveTo(sx(x),sy(y)); else ctx.lineTo(sx(x),sy(y));} ctx.stroke(); ctx.fillStyle='#1b2d55'; ctx.fillText('Smooth reference curve', W-180, 24); ctx.fillStyle='#d44500'; ctx.fillText('Line-segment fit', W-180, 44); document.getElementById('fitNote').innerText='Order '+order+' uses '+pieces+' line segments in this teaching display. Higher order follows local shape more closely, but it can also track noise.';")
statashiny, build export("webpage/nnsreg_interactive.html")

cd "webpage"
wdinit nnsreg_web, replace

webdoc put <h1>nnsreg: NNS regression in Stata</h1>
webdoc put <p><strong>nnsreg</strong> fits a univariate NNS-style piecewise linear curve and posts the result as a Stata estimator. Use it to inspect bends, flat regions, slope changes, and whether those features survive sensitivity checks. The examples here match the help file, README, and Stata Journal draft.</p>

webdoc put <h2>Questions It Helps Answer</h2>
webdoc put <p>Use <code>nnsreg</code> when you need to ask where a relationship bends, which part of the support is steep or flat, whether slope changes are stable across options, and whether a segmented fit adds value relative to <code>regress</code> or <code>npregress kernel</code>. Kernel smoothing remains a strong choice for smooth conditional-mean recovery; NNS is most useful when segmentation, clustered support, sharp local features, or stored slope-change output matter.</p>

webdoc put <h2>Syntax And Options</h2>
webdoc put <pre><code>nnsreg depvar indepvar [if] [in] [, order(integer) noise(mean|median) method(connect|ols) partition(xy|xonly) noplot graph_opts(string) generate(newvar)]</code></pre>
webdoc put <ul><li><code>order(#)</code> controls partition depth. Larger values allow more line segments.</li><li><code>noise(mean)</code> uses means for centers and regression points.</li><li><code>noise(median)</code> is a robustness check for skewness or outliers.</li><li><code>method(connect)</code> connects NNS regression points directly.</li><li><code>method(ols)</code> estimates an OLS spline at NNS-derived knots.</li><li><code>partition(xy)</code> is the default NNS-style quadrant partition.</li><li><code>partition(xonly)</code> is a vertical-band sensitivity check.</li></ul>

webdoc put <h2>Simulation Design</h2>
webdoc put <p>The education example simulates an S-shaped curve with diminishing returns. The water example simulates a flatter lower range and a sharper increase after an infrastructure-age threshold. The clustered benchmark simulates dense local regions with a peak, sharp dip, and recovery. These are simulated examples, so the figures test software behavior against known curves rather than estimating real policy effects.</p>

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
webdoc put <figure><img src="../manuscript/fig_edu_combined_compare.png" class="img-fluid" alt="Combined education comparison figure"><figcaption>Fit comparison and stored coefficients for the education simulation.</figcaption></figure>

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
webdoc put <figure><img src="../manuscript/fig_water_combined_compare.png" class="img-fluid" alt="Combined water comparison figure"><figcaption>Fit comparison and stored coefficients for the water simulation.</figcaption></figure>

webdoc put <h2>Worked Example: Clustered Benchmark</h2>
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
webdoc put <figure><img src="../manuscript/fig_cluster_combined_compare.png" class="img-fluid" alt="Combined clustered benchmark comparison figure"><figcaption>The clustered benchmark shows a case where NNS OLS improves recovery of a known curve relative to <code>npregress kernel</code>.</figcaption></figure>

webdoc put <h2>Diagnostics And Sensitivity</h2>
webdoc put <p>The verification script exports <code>table_model_diagnostics.csv</code> with a common in-sample <code>R^2</code>, observed RMSE, and true-curve RMSE. In the current run, <code>npregress kernel</code> and NNS OLS are close in the smooth education and water simulations, while NNS OLS improves substantially in the clustered peak/dip benchmark.</p>
webdoc put <figure><img src="../manuscript/fig_edu_sensitivity.png" class="img-fluid" alt="Education sensitivity figure"><figcaption>Education sensitivity to <code>order()</code> and <code>noise()</code>.</figcaption></figure>
webdoc put <figure><img src="../manuscript/fig_water_sensitivity.png" class="img-fluid" alt="Water sensitivity figure"><figcaption>Water sensitivity to <code>order()</code> and <code>noise()</code>.</figcaption></figure>
webdoc put <figure><img src="../manuscript/fig_cluster_sensitivity.png" class="img-fluid" alt="Clustered benchmark sensitivity figure"><figcaption>Clustered benchmark sensitivity to <code>order()</code>.</figcaption></figure>

webdoc put <h2>Interactive Fit Explorer</h2>
webdoc put <p>The statashiny companion below illustrates how increasing partition order can make a line-segment fit follow the local curve more closely while also increasing sensitivity to noise. It is a teaching visualization, not a separate Stata estimation result.</p>
webdoc put <div style="max-width:100%;overflow:hidden"><iframe src="nnsreg_interactive.html" style="width:100%;height:760px;border:1px solid #d8dee8;border-radius:8px;"></iframe></div>

webdoc put <h2>Install</h2>
webdoc put <pre><code>net install nnsreg, from("https://raw.githubusercontent.com/ericbooth/NNS-stata-public/main/github") replace</code></pre>

webdoc close
