# OBGlutamateDynamics
<img src="https://github.com/andsinpants/andsinpants.github.io/blob/main/static/images/papers/Moran2021a.PNG" alt="Moran2021a">
<img src="https://github.com/andsinpants/andsinpants.github.io/blob/main/static/images/papers/Moran2021b.PNG" alt="Moran2021b">


This is analysis code for Moran et al., 2021a and 2021b glutamate dynamics papers.
<ul><li><a href="https://www.eneuro.org/content/8/2/ENEURO.0110-21.2021" target="_blank">Moran et al., 2021a</a></li>
<li><a href="https://www.frontiersin.org/articles/10.3389/fncir.2021.779056/full" target="_blank">Moran et al., 2021b</a></li></ul>

 
<h2>The following <i><b><u>software</i></b></u> was utilized:</h2>
<dl><dt>MATLAB2018b</dt></dl>
<dd>Note: MATLAB2018b was the last edition used with this code, may not work on newer versions. </dd></dl>

<h2>The following <i><b><u>folders</i></b></u> are as follows:</h2>
<dl><dt>data</dt></dl>
<dd>Note: The data folders are not uploaded to github due to filesize. See below.</dd></dl>
<dl><dt>scripts</dt></dl>
<dl><dt>functions</dt></dl>

<h2>The following <i><b><u>functions</i></b></u> are as needed:</h2>
<dl><dt>autoArrangeFigures</dt>
<dd>This will autoarrange figures based on screen resolution for easy viewing (~<=27 figs @ 1920 x 1080).</dd></dl>
<dl><dt>cab</dt>
<dd>Close all but current figure dictated by figure handles.</dd></dl>
<dl><dt>dir2</dt>
<dd>Close all but current figure dictated by figure handles.</dd></dl>
<dl><dt>normalised_diff (custom)</dt>
<dd>Normalize between -1 and 1 based on the max ROI activity.</dd></dl>
<dl><dt>normrois (custom)</dt>
<dd>Normalize ROIs, requires normalised_diff function.</dd></dl>
<dl><dt>odorexp (custom)</dt>
<dd>To be removed. Ignore.</dd></dl>
<dl><dt>fitexp (custom)</dt>
<dd>Fit exponential to given points.</dd></dl>

<h2>Subfolders in data</h2>
<h3><b>Data was not uploaded to the github repo. <br><br>Ask <a href="mailto:andrew.moran@cuanschutz.com">Andrew Moran</a>
 or <a href="mailto:matt.wachowiak@utah.com">Matt Wachowiak</a>
 for data if you want to use for your own analysis.</b></h3>
<h3>2P</h3>
  <h4>odor_gun</h4>
    <ul>
    <li>pcd or tbet SF-iGluSnFR.A184S</li>
    <li>tbet SF-iGluSnFR.A184V awake</li>  
    <li>tbt iGluSnFR APV_NBQX only odor_gun</li>
    </ul>
  <h4>olfactometer</h4>
    <ul>
    <li>cck GCaMP CGP+APV_NBQX</li>
    <li>pcd GCaMP CGP + APV_NBQX</li>
    <li>pcd GCaMP CGP only</li>
    <li>pcd iGluSnFR 3x odors</li>
    <li>pcd or tbet SF-iGluSnFR.A184S</li>
    <li>pcd or tbet SF-iGluSnFR.A184V</li>
    <li>pcd or tbet SF-iGluSnFR.S72A</li>
    <li>tbet SF-iGluSnFR.A184V awake</li>
    <li>tbt iGluSnFR CGP + APV_NBQX APc-inj</li>
    </ul>
<h3>Epi</h3>
  <h4>olfactometer</h4>
    <ul><li>pcd iGluSnFR CGP and APV_NBQX</li></ul>
    
<h2>Subfolders in scripts</h2>
<h3>single color</h3>
  <ul>
  <li>pharmacology</li>
  <li>ITA</li>
  <li>odor_gun</li>
  </ul>
<h3>2-color</h3>
  <ul>
  <li>odor_gun</li>
  <li>ITA</li>
  </ul>