#! /usr/bin/python3
# Last edited on 2025-04-13 12:06:20 by stolfi
 
import os, sys, subprocess
from math import sqrt, hypot

pid = os.getpid()
tmp = f"/tmp/{pid}"

def main():
  testDoc = "cbabel WITHOUT multiscale"
  # run_test("cbabel", testDoc, (256,192), 1.000, 0.000, True, "zero", 0.8, 0)

  testDoc = "cbabel WITH multiscale"
  # run_test("cbabel", testDoc, (256,192), 1.000, 0.000, True, "zero", 0.8, 999)

  testDoc = "cplat3 WITHOUT synergy hintZ"
  # run_test("cplat3", testDoc, (256,192), 1.000, 0.000, True, "zero", 0.8, 999)

  testDoc = "cplat3 WITH synergy hintZ"
  # run_test("cplat3", testDoc, (256,192), 1.000, 0.005, True, "zero", 0.8, 999)

  testDoc = "fourg1 WITHOUT synergy hintZ"
  # run_test("fourg1", testDoc, (512,422), 1.000, 0.000, True, "zero", 0.8, 999)

  testDoc = "fourg1 WITH synergy hintZ"
  # run_test("fourg1", testDoc, (512,422), 1.000, 0.005, True, "zero", 0.8, 999)

  testDoc = "fours1 WITHOUT synergy hintZ"
  # run_test("fours1", testDoc, (512,384), 1.000, 0.000, True, "zero", 0.8, 999)

  testDoc = "fours1 WITH synergy hintZ"
  run_test("fours1", testDoc, (512,384), 1.000, 0.005, True, "zero", 0.8, 999)

  return 0
  # ......................................................................
  
def run_test(testName, testDoc, size, NGScale, hintZ_weight, refZ_use, initZ_source, initZ_noise, maxLevel):
  # Runs the recursive integration program on one set of images obtained by 
  # photometric (PM) stereo or multifocus (MF) stereo. 
  #
  # The {size} parameter should be a list {(NX,NY)} with the col and row counts 
  # of the input gradient or normal map.
  #
  # The {testName} is the name of the input dataset as in the {in} directory. The
  # {testDoc} is an explanatory string printed to {stderr}.
  #
  # The {NGScale} parameter is used when the input is a normal map. The
  # slopes will be {NGScale*ux/uz} and {NGScale*uy/uz} where
  # {(ux,uy,uz)} is the normal vector.  See the "-normals" option
  # of the integration program {prog}
  #
  # The {hintZ_weight} is the weight of the height hint (mutifocal Z)
  # map, relative to the gradient, in the synergistic integration. If it
  # is zero, no Z hint map is used, only the gradient. See the "-hint" option
  # of {prog}.
  #
  # The {refZ_use} specifies whether the input dataset includes a
  # reference heigh map for comparison with the computed height map. See
  # the "-refernce" option of {prog}.
  #
  # The {initZ_source} and {initZ_noise} options are the source of the
  # initial guess to be used ("zero", "hints", "reference") and the
  # relative amount of noise to add to it. See the "-initial" oprion of
  # {prog}.
  # 
  # The {maxLevel} integer is passed to {prog} as the "-maxLevel" option.
  #
  # Inputs, in directory {inDir}:
  #
  #   "PM-G.fni"  Gradient (slope) map from photometric stereo.
  #   "PM-N.fni"  Normal map from photometric stereo.
  #   "MF-Z.fni"  Height map from multifocus stereo.
  #   "RF-Z.fni"  Reference (ideal) height map for evaluation.
  #
  # The normal map is used only if the gradient map is not present.
  # The multifocus "MF-Z.fni" and reference RF-Z.fni height maps are optional.
  # The intergration will use "MF-Z.fni" as the height hint if the file exists 
  # and {hintZ_weight} is nonzero.
  #
  # Outputs, in folder {outDir}, for each two-digit {level} from 00 to the max level used:
  #
  #   "out-{level}-beg-G.fni"    Input gradient map (possibly from normal map), scaled.
  #   "out-{level}-beg-H.fni"    Input hints map, if given.
  #   "out-{level}-beg-R.fni"    Input reference map, if given.
  #   "out-{level}-beg-Z.fni"    Initial height map, with noise.
  #   "out-{level}-beg-E.fni"    Initial error map (beg-Z minus beg-R).
  #   "out-{level}-beg-E.txt"    One-line summary of initial error map.

  #   "out-{level}-beg-S.txt"    The level integation linear system derived from out-beg-G, out-beg-H.
  #   "out-{level}-beg-SW.fni"   Map with total weight of the equation for each pixel.

  #   "out-{level}-end-Z.fni"    Final height map.
  #   "out-{level}-end-E.fni"    Final error map (end-Z minus beg-R).
  #   "out-{level}-end-E.txt"    One-line summary of end-E.fni.
  
  # Also,in folder {outDir}/out-iters, for each two-digit {level} and some  
  # nine-digit iteration count {iter}:

  #   "out-iters/it-{level}-{iter}-Z.fni"    Height map after iteration {iter}.
  #   "out-iters/it-{level}-{iter}-E.fni"    Error map (Z.fni minus {outDir}/out-{level}-beg-R.fni).
  #   "out-iters/it-{level}-{iter}-E.txt"    One-line summary of E.fni.
  
  progDir = ".."
  prog = "gus_integrate_recursive"

  Er("starting {run_test} ...\n")
  Er(f"{testDoc}\n")

  xsize = f"{size[0]:04d}x{size[1]:04d}"  # Size for file names.
  inDir = f"in/{testName}/{xsize}"
  
  in_slopes_fni = f"{inDir}/PM-G.fni"
  in_normals_fni = f"{inDir}/PM-N.fni"
  in_hintZ_fni = f"{inDir}/MF-Z.fni"
  in_refZ_fni = f"{inDir}/RF-Z.fni"
    
  # Shall it use a hint height map:
  hintZ_use = file_OK(in_hintZ_fni) and hintZ_weight > 0

  xhintZ = "hintT" if hintZ_use else "hintF" # {hintZ_use} for file names."
  outDir = f"out/{testName}/{xsize}-{xhintZ}"
  bash(f"mkdir -p {outDir}")
  outPrefix = f"{outDir}/out"
  
  iterDir = f"{outPrefix}-iters"
  iterPrefix = f"{iterDir}/it"

  # Input slope or normal map otpion:
  if file_OK(in_slopes_fni):
    input_options = [ "-slopes", in_slopes_fni, ]
  elif file_OK(in_normals_fni):
    input_options = [ "-normals", in_normals_fni, "scale", f"{NGScale:6f}", ]
  else:
    assert False, f"{in_slopes_fni} and {in_normals_fni} are both missing"

  # Reference height map otpion for error monitoring:
  if refZ_use:
    out_errZ_fni = f"{outPrefix}-00-end-E.fni"
    refZ_options = [ "-reference", in_refZ_fni, ]
  else:
    out_errZ_fni = None
    refZ_options = [ ]

  # Hint height map option:
  if hintZ_use:
    hintZ_options = [ "-hints", in_hintZ_fni, f"{hintZ_weight:.6f}", ]
  else:
    hintZ_options = [ ]
    
  Er("showing input files ...\n")
  for state, tag, fniFile, viewScale in \
      ( ( 'input',         'PM-G', in_slopes_fni, "1.000" ),
        ( 'input',         'PM-N', in_normals_fni, "auto" ),
        ( 'input',         'MF-Z', in_hintZ_fni, "1.000" ),
        ( 'input',         'RF-Z', in_refZ_fni, "1.000" ),
      ):
    if file_OK(fniFile):
      showfni('input', tag, fniFile, viewScale)

  initZ_options = [ "-initial", initZ_source, f"{initZ_noise:6f}", ]

  # Clear the gray level scale:
  if size == (512, 422):
    clear_options = [ "-clear", "0", "511", "0", "40", ]
  else:
    clear_options = [ ]

  bash("rm -f {%s,%s}*.{fni,sys,pgm,txt,png}" % (outPrefix, iterPrefix))
  cmd = \
    [ f"{progDir}/{prog}" ] \
    + input_options \
    + hintZ_options  \
    + initZ_options \
    + refZ_options \
    + clear_options \
    + [ "-maxLevel", f"{maxLevel:d}",
        "-convTol", "0.0005",
        "-maxIter", "250",
        "-outPrefix", outPrefix,
        "-reportStep", "10",
        "-verbose"
      ]
  Er("running command\n")
  Er(" ".join(cmd) + "\n")
  
  run_command(cmd);

  Er("showing results ...\n")
  # Initial guess at level 0:
  out_iniZ_fni = f"{outPrefix}-ini-Z.fni"
  showfni( 'initial guess', 'Z', out_iniZ_fni, "auto" )
  for kl in range(20):
    level = 20 - kl
    levelPrefix = f"{outPrefix}-{level:02d}"
    for state in ("beg", "end" ):
      for tag in ( 'Z', 'G', 'E' ):
        if tag == "E":
          viewScale = "auto";
        else:
          viewScale = "1.000"
        fniFile = f"{levelPrefix}-{state}-{tag}.fni"
        Er(f"showing {fniFile} ...\n")
        showfni( state, tag, fniFile, viewScale )

  make_movies(outDir, iterDir, size, hintZ_use)
  
  Er("done {run_test}.\n")
  return
  # ......................................................................
  
def make_movies(outDir, iterDir, size, hintZ_use):  
  
  Er("making the heigh map animation ...\n")
  zMax = hypot(size[0], size[1])/sqrt(8)
  bash(f"make_iteration_movie.sh {outDir} {iterDir} 10 Z -{zMax:.4f} +{zMax:.4f}")
  
  Er("making the height error map animation ... \n");
  eMax = 2.0 if hintZ_use else zMax
  bash(f"make_iteration_movie.sh {outDir} {iterDir} 10 E -{eMax:.8f} +{eMax:.8f}")
  
  return
  # ......................................................................

def showfni(state, tag, fniFile, viewScale):
  fniName = os.path.basename(fniFile)

  stg = f"{state} {tag}"

  if file_OK(fniFile):
    Er(f"showing {fniFile} ...\n")

    Er(f"state = '{state}'  tag = '{tag}'  scale = {viewScale}\n")
    title = f"{state} {tag} : {fniName}"
    pgm_file = f"{tmp}.pgm"
    png_file = "${fniFile/.fni/.png}"
    bash(f"fni_view -title '{title}' -scale {viewScale} {fniFile}")
    # bash(f"fni_to_pnm -channel 0 -yAxis up < {fniFile} > {pgm_file}")
    # bash(f"convert {pgm_file} {png_file}")
    # bash(f"display -title '{title}' -filter box -resize 'x800<' {png_file}")
  else:
    Er(f"!! {fniFile} not present\n")
  return
  # ......................................................................

def Er(msg):
  sys.stderr.write(msg)
  # ......................................................................
  
def bash(cmd):
  # The {cmd} must be a string with a command in {bash} syntax. 
  # Executes {cmd} with "/bin/bash".
  result = subprocess.run([ cmd ], shell=True, executable="/bin/bash")
  print(result.stderr)
  print(result.stdout)
  assert result.returncode == 0, f"** {cmd[0]} failed - returned status = {result}"
  return
  # ......................................................................

def file_OK(path):
  return path != None and os.path.exists(path) and os.path.getsize(path) > 0
  # ......................................................................

def run_command(cmd):
  # The {cmd} must be a list where the first elem is an executable file name
  # and the other elements are the command line arguments to be passed 
  # to its as the {argv} array. 
  result = subprocess.run(cmd, text = True)
  print(result.stderr)
  print(result.stdout)
  assert result.returncode == 0, f"** {cmd[0]} failed - returned status = {result}"
  return
  # ......................................................................
  
main()
