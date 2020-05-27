function select_fct() {
  var val = document.getElementById("access").selectedIndex;
  var x = document.getElementById("devBlock");
  if (val == 1) {
    x.style.display = "block";
  } else {
    x.style.display = "none";
  }
  generate_all();
}

function generate_all() {
  var access = document.getElementById("access").value;
  var cuda_path = document.getElementById("cuda_path").value;
  var conda_path = document.getElementById("conda_path").value;
  var magma_install = document.getElementById("magma_install").value;
  var cuda_sm = document.getElementById("cuda_sm").value;
  var hf16 = "OFF";
  if (document.getElementById("half16").checked) {
    hf16 = "ON";
  }

  generate_bashrc(
    document.getElementById("bashrc"),
    access,
    conda_path,
    cuda_path,
    magma_install,
    hf16
  );
  generate_script(document.getElementById("script"), access, cuda_sm);
}

function generate_bashrc(
  div_text,
  access,
  conda_path,
  cuda_path,
  magma_install,
  hf16
) {
  text_html = `
  <p><h3>Add to your .bashrc</h3></p>
    <div class="fragment">

  <div class="line">## CONDA default definitions</div>
  <div class="line">export CONDA_ROOT=${conda_path}</div>
  <div class="line">export PATH=$CONDA_ROOT/bin:$PATH</div>
  <div class="line"></div>
  `;
  if (access == "conda") {
    text_html += `
    <div class="line">#COMPASS default definitions</div>
    <div class="line">export SHESHA_ROOT=$HOME/shesha</div>
    <div class="line">export PYTHONPATH=$NAGA_ROOT:$SHESHA_ROOT:$PYTHONPATH</div>
    `;
  } else {
    text_html += `
    <div class="line">## CUDA default definitions</div>
    <div class="line">export CUDA_ROOT=${cuda_path}</div>
    <div class="line">export CUDA_INC_PATH=$CUDA_ROOT/include</div>
    <div class="line">export CUDA_LIB_PATH=$CUDA_ROOT/lib</div>
    <div class="line">export CUDA_LIB_PATH_64=$CUDA_ROOT/lib64</div>
    <div class="line">export PATH=$CUDA_ROOT/bin:$PATH</div>
    <div class="line">export LD_LIBRARY_PATH=$CUDA_LIB_PATH_64:$CUDA_LIB_PATH:$LD_LIBRARY_PATH</div>
    <div class="line"></div>
    <div class="line">#MAGMA definitions</div>
      <div class="line">export MAGMA_ROOT="${magma_install}"</div>
      <div class="line">export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MAGMA_ROOT/lib</div>
      <div class="line">export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$MAGMA_ROOT/lib/pkgconfig</div>
    <div class="line"></div>
      <div class="line">#COMPASS default definitions</div>
      <div class="line">export COMPASS_ROOT=$HOME/compass</div>
      <div class="line">export COMPASS_INSTALL_ROOT=$COMPASS_ROOT/local</div>
      <div class="line">export COMPASS_DO_HALF="${hf16}" # set to ON if you want to use half precision RTC (needs SM>=60)</div>
      <div class="line">export NAGA_ROOT=$COMPASS_ROOT/naga</div>
      <div class="line">export SHESHA_ROOT=$COMPASS_ROOT/shesha</div>
      <div class="line">export LD_LIBRARY_PATH=$COMPASS_INSTALL_ROOT/lib:$LD_LIBRARY_PATH</div>
      <div class="line">export PYTHONPATH=$NAGA_ROOT:$SHESHA_ROOT:$COMPASS_INSTALL_ROOT/python:$PYTHONPATH</div>
      <div class="line">export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$COMPASS_INSTALL_ROOT/lib/pkgconfig</div>
    <div class="line"></div>
    <div class="line">#third party lib path</div>
      <div class="line">export CUB_ROOT=$COMPASS_ROOT/tplib/cub</div>
      <div class="line">export WYRM_ROOT=$COMPASS_ROOT/tplib/wyrm</div>
    `;
  }
  text_html += "</div><!-- fragment -->";

  div_text.innerHTML = text_html;
}

function generate_script(div_text, access, cuda_sm) {
  // var val = document.getElementById("access").selectedIndex;
  // var x = document.getElementById("devBlock");
  // if (val == 1) {
  //   x.style.display = "block";
  // } else {
  //   x.style.display = "none";
  // }
  text_html = `
  <p><h3>Script to run</h3></p>
  <div class="fragment">
  <div class="line">source $HOME/.bashrc</div>
  <div class="line"></div>
  <div class="line">mkdir -p $HOME/tmp_compass</div>
  <div class="line">cd $HOME/tmp_compass</div>
  <div class="line"></div>
  <div class="line">wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh</div>
  <div class="line">bash Miniconda3-latest-Linux-x86_64.sh -b -p $CONDA_ROOT</div>
  <div class="line"></div>
  `;
  if (access == "conda") {
    text_html += `
    <div class="line">conda install -y -c compass compass</div>
      <div class="line">cd $HOME</div>
      <div class="line">git clone https://github.com/ANR-COMPASS/shesha.git</div>
    `;
  } else {
    text_html += `
      <div class="line">export MKLROOT=$CONDA_ROOT</div>
      <div class="line">export CUDADIR=$CUDA_ROOT</div>
      <div class="line">export NCPUS=8</div>
      <div class="line">export GPU_TARGET=sm_${cuda_sm}</div>
      <div class="line"></div>
      <div class="line">conda install -y numpy mkl-devel pyqtgraph ipython pyqt qt matplotlib astropy blaze h5py hdf5 pytest pandas scipy docopt tqdm tabulate</div>
      <div class="line"></div>
      <div class="line">wget http://icl.cs.utk.edu/projectsfiles/magma/downloads/magma-2.5.2.tar.gz -O - | tar xz</div>
      <div class="line">cd magma-2.5.2</div>
      <div class="line"></div>
      <div class="line">cp make.inc-examples/make.inc.mkl-gcc make.inc</div>
      <div class="line">sed -i -e "s:/intel64: -Wl,-rpath=$CUDADIR/lib64 -Wl,-rpath=$MKLROOT/lib:" make.inc</div>
      <div class="line"></div>
      <div class="line">make -j $NCPUS shared sparse-shared</div>
      <div class="line">make install prefix=$MAGMA_ROOT</div>
      <div class="line"></div>
      <div class="line">cd $HOME</div>
      <div class="line">git clone https://gitlab.obspm.fr/compass/compass --recurse-submodules</div>
      <div class="line">cd $COMPASS_ROOT</div>
      <div class="line">./compile.sh</div>
    `;
  }
  text_html += "</div><!-- fragment -->";

  div_text.innerHTML = text_html;
}

function generate_form(div_text) {
  div_text.innerHTML = `
  <fieldset>
    <legend>Option script</legend>\n
    <label for="access">Access : </label>\n
    <select name="access" id="access" onchange="select_fct()">\n
      <option value="conda" selected="selected">Conda</option>\n
      <option value="git">Gitlab Obspm</option>\n
    </select><br>\n
    <label for="conda_path">Conda install:</label>\n
    <input type="text" name="conda_path" id="conda_path" value="$HOME/miniconda3" onchange="generate_all()" /><br>
    <div id="devBlock">
      <label for="cuda_path">Cuda install:</label>
      <input type="text" name="cuda_path" id="cuda_path" value="/usr/local/cuda" onchange="generate_all()" /><br>
      <label for="magma_install">Magma installation path:</label>
      <input type="text" name="magma_install" id="magma_install" value="$HOME/local/magma" onchange="generate_all()" /><br>
      <label for="cuda_sm">CUDA Compute capability:</label>
      <input type="text" name="cuda_sm" id="cuda_sm" value="70" onchange="generate_all()" /><br>
      <label for="half16">half16</label>
      <input type="checkbox" name="half16" id="half16" onchange="generate_all()" /><br>
    </div>
  </fieldset>
  <div id="bashrc">
  </div>
  <div id="script">
  </div>
  `;
}

generate_form(document.getElementsByClassName("textblock")[0]);
document.getElementById("devBlock").style.display = "none";
generate_all();
