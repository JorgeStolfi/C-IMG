#! /usr/bin/python3
# Last edited on 2025-04-12 10:40:47 by stolfi
 
import os, sys, subprocess
from math import sqrt, hypot

pid = os.getpid()
tmp = f"/tmp/{pid}"

def Er(msg):
  sys.stderr.write(msg)
  # ......................................................................
  
def bash(cmd):
  # The {cmd} must be a string with a command in {bash} syntax. 
  # Executes {cmd} with "/bin/bash".
  subprocess.run([ cmd ], shell=True, executable="/bin/bash")
  # ......................................................................

def file_ok(path):
  return os.path.exists(path) and os.path.getsize(path) > 0
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

def main():
  testDoc = "cbabel WITHOUT multiscale"
  # run_test("cbabel", testDoc, (256,192), 1.000, 1.000, 0.000, True, "ZR", 0.8, 0)

  testDoc = "cbabel WITH multiscale"
  # run_test("cbabel", testDoc, (256,192), 1.000, 1.000, 0.000, True, "ZR", 0.8, 999)

  testDoc = "cplat3 WITHOUT synergy hint"
  # run_test("cplat3", testDoc, (256,192), 1.000, 1.000, 0.000, True, "ZR", 0.8, 999)

  testDoc = "cplat3 WITH synergy hint"
  # run_test("cplat3", testDoc, (256,192), 1.000, 1.000, 0.005, True, "ZR", 0.8, 999)

  testDoc = "fourg1 WITHOUT synergy hint"
  # run_test("fourg1", testDoc, (512,422), 1.000, 1.000, 0.000, True, "ZR", 0.8, 999)

  testDoc = "fourg1 WITH synergy hint"
  # run_test("fourg1", testDoc, (512,422), 1.000, 1.000, 0.005, True, "ZR", 0.8, 999)

  testDoc = "fours1 WITHOUT synergy hint"
  # run_test("fours1", testDoc, (512,422), 1.000, 1.000, 0.000, True, "ZR", 0.8, 999)

  testDoc = "fours1 WITH synergy hint"
  run_test("fours1", testDoc, (512,384), 1.000, 1.000, 0.005, True, "ZR", 0.8, 999)

  return 0
  # ......................................................................
  
def run_test(testName, testDoc, size, GScale, ZScale, hintsWeight, hasRefZ, initial_heights, initial_noise, maxLevel):
  # Runs the recursive integration program on one set of images obtained by 
  # photometric (PM) stereo or multifocus (MF) stereo. 
  #
  # Inputs, in directory {inDir}:
  #
  #   PM-G.fni  Gradient (slope) map from photometric stereo.
  #   PM-N.fni  Normal map from photometric stereo.
  #   MF-Z.fni  Height map from multifocus stereo.
  #   RF-Z.fni  Reference (ideal) height map for evaluation.
  #
  # The normal map is used only if the gradient map is not present.
  # The multifocus and reference height maps are optional.
  #
  # Outputs, in folder {outDir}, for each two-digit {level} from 00 to the max level used:
  #
  #   out-{level}-beg-G.fni    Input gradient map (possibly from normal map), scaled.
  #   out-{level}-beg-H.fni    Input hints map, if given.
  #   out-{level}-beg-R.fni    Input reference map, if given.
  #   out-{level}-beg-Z.fni    Initial height map, with noise.
  #   out-{level}-beg-E.fni    Initial error map (beg-Z minus beg-R).
  #   out-{level}-beg-E.txt    One-line summary of initial error map.

  #   out-{level}-beg-S.txt    The level integation linear system derived from out-beg-G, out-beg-H.
  #   out-{level}-beg-SW.fni   Map with total weight of the equation for each pixel.

  #   out-{level}-end-Z.fni    Final height map.
  #   out-{level}-end-E.fni    Final error map (end-Z minus beg-R).
  #   out-{level}-end-E.txt    One-line summary of end-E.fni.
  
  # Also,in folder {outDir}/out-iters, for each two-digit {level} and some  
  # nine-digit iteration count {iter}:

  #   out-iters/it-{level}-{iter}-Z.fni    Height map after iteration {iter}.
  #   out-iters/it-{level}-{iter}-E.fni    Error map (Z.fni minus {outDir}/out-{level}-beg-R.fni).
  #   out-iters/it-{level}-{iter}-E.txt    One-line summary of E.fni.
  
  progDir = ".."
  prog = "gus_integrate_recursive"

  Er("starting {run_test} ...\n")
  Er(f"{testDoc}\n")

  xsize = f"{size[0]:04d}x{size[1]:04d}"
  inDir = f"in/{testName}/{xsize}"

  outDir = f"out/{testName}/{xsize}"
  bash(f"mkdir -p {outDir}")
  outPrefix = f"{outDir}/out"
  
  iterDir = f"{outPrefix}-iters"
  iterPrefix = f"{iterDir}/it"

  in_slopes_fni = f"{inDir}/PM-G.fni"
  in_normals_fni = f"{inDir}/PM-N.fni"
  if file_ok(in_slopes_fni):
    input_options = [ "-slopes", in_slopes_fni, "scale", f"{GScale:6f}", ]
  elif file_ok(in_normals_fni):
    input_options = [ "-normals", in_normals_fni, ]
  else:
    assert False, f"{in_slopes_fni} and {in_normals_fni} are both missing"

  if hasRefZ:
    in_refz_fni = f"{inDir}/RF-Z.fni"
    out_errz_fni = f"{outPrefix}-00-end-E.fni"
    reference_options = [ "-reference", in_refz_fni, ]
  else:
    in_refz_fni = None
    out_errz_fni = None
    reference_options = [ ]

  if hintsWeight != 0:
    in_hints_fni = f"{inDir}/MF-Z.fni"
    hints_options = [ "-hints", in_hints_fni, f"{hintsWeight:.6f}", ]
  else:
    in_hints_fni = None
    hints_options = [ ]

  Er("showing input files ...\n")
  viewZScale = ZScale*2
  for tag, fniFile, viewScale in \
      ( ( 'PM-G', in_slopes_fni, f"{GScale:.4f}" ),
        ( 'PM-N', in_normals_fni, "auto" ),
        ( 'MF-Z', in_hints_fni, f"{viewZScale:.4f}" ),
        ( 'RF-Z', in_refz_fni, f"{viewZScale:.4f}" )
      ):
    showfni('input', tag, fniFile, viewScale)

  if initial_heights == "ZR":
    initial_options = [ "-initial", "zero", f"{initial_noise:6f}", ]
  elif initial_heights  == "MF" :
    initial_options = [ "-initial", "hints", f"{initial_noise:6f}", ]
  elif initial_heights  == "RF" :
    initial_options = [ "-initial", "reference", f"{initial_noise:6f}", ]
  else:
    assert False, f"invalid initial_heights = '{initial_heights}'"

  # Clear the gray level scale:
  if size == (512, 422):
    clear_options = [ "-clear", "0", "511", "0", "40", ]
  else:
    clear_options = [ ]

  bash("rm -f {%s,%s}*.{fni,sys,pgm,txt,png}" % (outPrefix, iterPrefix))
  cmd = \
    [ f"{progDir}/{prog}" ] \
    + input_options \
    + hints_options  \
    + initial_options \
    + reference_options \
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
  for kl in range(20):
    level = 20 - kl
    levelPrefix = f"{iterPrefix}-{level:02d}"
    for state in ("beg", "end" ):
      for tag in ( 'Z', 'G', 'E' ):
        if tag == "G":
          viewScale = "auto";
        else:
          viewScale = "1.000"
        fniFile = f"{levelPrefix}-{state}-{tag}.fni"
        showfni( state, tag, fniFile, viewScale )

  make_movies(outDir, iterDir, size)
  
  Er("done {run_test}.\n")
  return
  # ......................................................................
  
def make_movies(outDir, iterDir, size):  
  
  Er("making the heigh map animation ...\n")
  zMax = hypot(size[0], size[1])/sqrt(8)
  bash(f"make_iteration_movie.sh {outDir} {iterDir} 10 Z -{zMax:.4f} +{zMax:.4f}")
  
  Er("making the height error map animation ... \n");
  eMax = 2.0
  bash(f"make_iteration_movie.sh {outDir} {iterDir} 10 E -{eMax:.8f} +{eMax:.8f}")
  
  return
  # ......................................................................

def showfni(state, tag, fniFile, viewScale):
  fniName = os.path.basename(fniFile)

  stg = f"{state} {tag}"

  if file_ok(fniFile):
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
  
main()
