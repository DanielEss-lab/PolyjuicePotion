<h3>Polyjuice</h3>\
Created by Zack Meyer, Nathan Morgan, and Billy Hirschi

<p>Usage: Identifies monodentate ligands in an xyz file containing an organometallic
molecule with a single metal core and replaces those ligands with a methyl group.
This program can be used to replace a single monodentate ligand or all monodentate
ligands in the compound, however changing it does require a small edit to the source code
The source code also contains code for saving the compounds during various stages of the
replacement/identification process, all labelled with comments describing what each section
of code does.</p>

<p>Requirements: requires Python 3.6, openbabel 3.1.1, all .py files in the same directory,
the xyz files with the organometallic compound(s) in the same directory as the .py files<br>
Recommended: Anaconda</p>
<p>Install openbabel:<br><code>conda create --name myenv python=3.6</code> with
<code>myenv</code> being the desired environment name<br>
<code>conda activate myenv</code><br>
<code>conda install -c conda-forge openbabel</code></p>

<h4>IMPORTANT</h4>
<p>Make sure to run Polyjuice in the correct anaconda environment.</p>
<code>python3 main.py</code> to run Polyjuice.
    