function select_fct() {
    var val = document.getElementById('acces').selectedIndex;
    var x = document.getElementById('devBlock');
    if (val == 1) {
      x.style.display = 'block';
    } else {
      x.style.display = 'none';
    }
  }

document.write('<div class="textblock">');
document.write('  <form method="post" action="../../generate_script.php">');
document.write('    <fieldset>');
document.write('      <legend>Option script</legend>');
document.write('      <label for="acces">Accès : </label>');
document.write('      <select name="acces" id="acces" onchange="select_fct()">');
document.write('        <option value="conda" selected="selected">Conda</option>');
document.write('        <option value="git">Gitlab Obspm</option>');
document.write('      </select><br>');
document.write('      <label for="conda_path">Conda install:</label>');
document.write('      <input type="text" name="conda_path" id="conda_path" value="$HOME/miniconda3" /><br>');
document.write('      <div id="devBlock">');
document.write('        <label for="cuda_path">Cuda install:</label>');
document.write('        <input type="text" name="cuda_path" id="cuda_path" value="/usr/local/cuda" /><br>');
document.write('        <label for="magma_install">Magma installation path:</label>');
document.write('        <input type="text" name="magma_install" id="magma_install" value="$HOME/local/magma" /><br>');
document.write('        <label for="cuda_sm">CUDA Compute capability:</label>');
document.write('        <input type="text" name="cuda_sm" id="cuda_sm" value="70" /><br>');
document.write('        <input type="checkbox" name="half16" id="half16" /> <label for="half16">half16</label><br>');
document.write('      </div>');
document.write('    </fieldset>');
document.write('    <input type="submit" value="envoyer" />');
document.write('  </form>');
document.write('</div>');

document.getElementById('devBlock').style.display = 'none';
