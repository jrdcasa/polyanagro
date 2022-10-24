import os
import sys
import re
import shutil
import glob
import argparse
from collections import defaultdict
import utils
import datetime
import numpy as np

# =============================================================================
def parse_arguments():

    def str2bool(v):
        if isinstance(v, bool):
            return v
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')



    desc = """Simple analysis of the VOTCA IBI results.
    This is part of the polyanagro library"""

    parser = argparse.ArgumentParser(description=desc)
    subparser = parser.add_subparsers(dest='command', required=True)

    ibi = subparser.add_parser('ibi')

    ibi.add_argument("-p", "--path", dest="path",
                        help="Path to the directory containing the step_* folders.",
                        action="store", required=True, default=None)

    ibi.add_argument("-b", "--begin", dest="begin",
                        help="Begin step. By default 0",
                        action="store", required=False, default=None)

    ibi.add_argument("-t", "--temp", dest="temp",
                        help="Temperature in K. By default 500K.",
                        action="store", required=False, default=500)

    ibi.add_argument("-e", "--end", dest="end",
                        help="End step. By default the last one presents in the path folder",
                        action="store", required=False, default=None)

    ibi.add_argument("--tmpdir", dest="tmpdir",
                        help="Temporal directory where png of each step is stored",
                        action="store_true", required=False, default=False)

    ibi.add_argument("--log", dest="log",
                        help="Name of the file to write logs from this command",
                        action="store", required=False, default="votca_analysis.log")

    ibi .add_argument("--press", dest="press", nargs=1,
                        help="Extract pressure from edr files in each step. "
                             "The argument for these argument is the path to gmx command of GROMACS",
                        action="store", required=False, default=None)

    args = parser.parse_args()

    subdir = os.path.join(args.path, "step_*")
    dirs = sorted(glob.glob(subdir))

    if len(dirs) == 0:
        print("Error!!! step_* folders from VOTCA calculation cannot be found!!.")
        print("Error!!! path: {} either does not exist or not a VOTCA folder!!".format(args.path))
        exit()

    if args.begin is None:
        args.begin = int(dirs[0].split("_")[-1])

    if args.end is None:
        args.end = int(dirs[-1].split("_")[-1])

    # Check that GROMACS gmx tool can be found.
    if args.press is not None:
        isgmx = os.path.isfile(args.press)
        if not isgmx:
            print("Error!!! Path to gmx does not exist in the system.")
            print("Error!!! path: {}".format(args.press))
            exit()

    kb_kjmol = 8.314462618e-3    # kJ/mol
    args.kbT = float(args.temp) * kb_kjmol

    return args

# =============================================================================
def print_header(version, logger_log=None):

    msg = """
    ***********************************************************************
                           Votca Analysis Tool
              ----------------------------------------------

                                Version {}

                              Dr. Javier Ramos
                      Macromolecular Physics Department
                Instituto de Estructura de la Materia (IEM-CSIC)
                               Madrid (Spain)

        This utility is part of the polyanagro library. Polyanagro is an 
        open-source python library to analyze simulations of polymer systems.

        This software is distributed under the terms of the
        GNU General Public License v3.0 (GNU GPLv3). A copy of
        the license (LICENSE.txt) is included with this distribution.

    ***********************************************************************
        """.format(version)

    print(msg) if logger_log is None else logger_log.info(msg)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

    m = "\t\tStart Job at {} ============".format(now)
    print(m) if logger_log is None else logger_log.info(m)

    m1 = ""
    for item in sys.argv[1:]:
        m1 += " {}".format(item)
    m = "\n\t\tCommand line: \n"
    m += "\t\t\tpython {}".format(os.path.split(sys.argv[0])[1])
    m += m1 + "\n"
    m += "\t\t\t         or\n"
    m += "\t\t\tvotca_analysis".format(os.path.split(sys.argv[0])[1])
    m += m1 + "\n"
    print(m) if logger_log is None else logger_log.info(m)

# =============================================================================
def try_to_guess_interaction(string):

    """
    Guess the type interactions by counting the number of "-" in the string
    Returns:

    """
    n = string.count("-")
    if n == 1:
        return "bond"
    elif n == 2:
        return "bend"
    elif n == 3:
        return "dihedral"
    else:
        return "non-bonded"

# =============================================================================
def read_distribution_files(args):

    # Target distributions from the first step
    pattern = os.path.join(args.path, "step_000/*.tgt".format(args.begin))
    tgt_list = sorted(glob.glob(pattern))
    pattern = os.path.join(args.path, "step_000/*.pot.new".format(args.begin))
    pot_list = sorted(glob.glob(pattern))

    if len(tgt_list) < 1:
        return None, None

    # Read target distributions
    target_distributions = defaultdict(dict)
    type_distributions_dict = defaultdict(list)
    for idx, itgt in enumerate(tgt_list):
        itgt_base = os.path.basename(itgt).split(".")[0]
        ipot = pot_list[idx]
        target_distributions[itgt_base] = {"xdata": [], "ydata": [], "type": "",
                                           "title": "", "filename":itgt, 'potfilename':ipot}
        title = os.path.basename(itgt).split(".")[0]
        typeinter = try_to_guess_interaction(title)
        type_distributions_dict[typeinter].append(title)
        with open(itgt, 'r') as f:
            lines = f.readlines()
            for iline in lines:
                if re.match("^#", iline):
                    continue
                x, y, _ = iline.split()
                target_distributions[itgt_base]["xdata"].append(float(x))
                target_distributions[itgt_base]["ydata"].append(float(y))
            target_distributions[itgt_base]["type"] = typeinter
            target_distributions[itgt_base]["title"] = title

    return target_distributions, type_distributions_dict

# =============================================================================
def prepare_individualdist_gnuplots(target_distributions, path_steps, terminal="png"):

    """
    Prepare gnu and png files for all distributions each step
    The png figures are located in the "path_steps+png" directory.

    Args:
        target_distributions: Dictionary with the info read from the VOTCA simulation
        path_steps: path where the steps_??? directories are located
        terminal: terminal line for gnuplot

    """

    if terminal == "png":
        terminal = "set terminal pngcairo size 550,400 enhanced font 'Verdana,10'\n" \
                   "set encoding iso_8859_1"
    else:
        terminal = terminal

    for itarget, ivalue in target_distributions.items():
        filename = ivalue["filename"]
        type_int = ivalue["type"]
        xmax = max(target_distributions[itarget]["xdata"])
        xmin = min(target_distributions[itarget]["xdata"])
        deltax = target_distributions[itarget]["xdata"][2] - target_distributions[itarget]["xdata"][1]
        ymax = max(target_distributions[itarget]["ydata"])
        ymin = min(target_distributions[itarget]["ydata"])
        deltay = (ymax - ymin)/10

        label = os.path.basename(filename).split(".")[0]
        fnamegnu = "ibi_dist_"+label+".gnu"
        pattern_step = os.path.join(path_steps, "step_*")

        with open(fnamegnu, 'w') as fgnu:
            fgnu.writelines("reset\n")
            fgnu.writelines(terminal)
            fgnu.writelines("\n")
            fgnu.writelines("#Target distribution\n")
            fgnu.writelines('TGT="{}"\n'.format(filename))
            fgnu.writelines("\n")
            fgnu.writelines('DIR=system("ls -d {0:s} |egrep -v step_000")\n'.format(pattern_step))
            fgnu.writelines("\n")
            fgnu.writelines('set title "{} beads"\n'.format(label))
            if type_int == "bond":
                fgnu.writelines('set xlabel "r (nm)"\n')
                fgnu.writelines('set ylabel "P_{bond}(nm{^-1})"\n')
                fgnu.writelines('set format x "%.1f"\n')
                fgnu.writelines('set format y "%.1f"\n')
                fgnu.writelines('set xrange [{}:{}]\n'.format(xmin-2*deltax, xmax+2*deltax))
                fgnu.writelines('set yrange [{}:{}]\n'.format(ymin, ymax+deltay))
                fgnu.writelines('set grid\n')
            elif type_int == "bend":
                fgnu.writelines('set xlabel "{/Symbol q} (rad)"\n')
                fgnu.writelines('set ylabel "P_{angle}(rad{^-1})"\n')
                fgnu.writelines('set format x "%.1f"\n')
                fgnu.writelines('set format y "%.1f"\n')
                fgnu.writelines('set xrange [{}:{}]\n'.format(xmin, xmax))
                fgnu.writelines('set yrange [{}:{}]\n'.format(ymin, ymax+deltay))
                fgnu.writelines('set grid\n')
            elif type_int == "dihedral":
                fgnu.writelines('set xlabel "{/Symbol f} (rad)"\n')
                fgnu.writelines('set ylabel "P_{dih}(rad{^-1})"\n')
                fgnu.writelines('set format x "%.3f"\n')
                fgnu.writelines('set format y "%.3f"\n')
                fgnu.writelines('set xrange [{}:{}]\n'.format(xmin, xmax))
                fgnu.writelines('set yrange [{}:{}]\n'.format(ymin, ymax+deltay))
                fgnu.writelines('set grid\n')
            elif type_int == "non-bonded":
                fgnu.writelines('set xlabel "r (nm)"\n')
                fgnu.writelines('set ylabel "g_{{{0:s}}}(r)"\n'.format(label))
                fgnu.writelines('set format x "%.1f"\n')
                fgnu.writelines('set format y "%.1f"\n')
                fgnu.writelines('set xrange [{}:{}]\n'.format(xmin, xmax))
                fgnu.writelines('set yrange [{}:{}]\n'.format(ymin, ymax+deltay))
                fgnu.writelines('set grid\n')
            fgnu.writelines('n=0\n')
            fgnu.writelines('do for [data in DIR] {\n')
            fgnu.writelines('   n=n+1\n')
            fgnu.writelines('   f1 = data."/{0:s}.dist.new"\n'.format(label))
            if type_int == "bond":
                fgnu.writelines('   set output sprintf("./png_dist/bond_{}_ibi%03d.png",n)\n'.format(label))
            if type_int == "bend":
                fgnu.writelines('   set output sprintf("./png_dist/bend_{}_ibi%03d.png",n)\n'.format(label))
            if type_int == "dihedral":
                fgnu.writelines('   set output sprintf("./png_dist/dihedral_{}_ibi%03d.png",n)\n'.format(label))
            if type_int == "non-bonded":
                fgnu.writelines('   set output sprintf("./png_dist/NB_{}_ibi%03d.png",n)\n'.format(label))
            if type_int == "bond" or type_int == "bend" or type_int == "dihedral":
                fgnu.writelines('   #conv = system("cat ".data."/{0:s}.conv")\n'.format(label))
                fgnu.writelines('   #merit = system("head -".n." merit.dat|tail -1")\n'.format(label))
                fgnu.writelines('   #set label "conv=".conv at 1.3,0.6\n')
                fgnu.writelines('   #set label "merit=".merit at 1.3,0.4\n')
            else:
                fgnu.writelines('   conv = system("cat ".data."/{0:s}.conv")\n'.format(label))
                fgnu.writelines('   merit = system("head -".n." merit_{}.dat|tail -1")\n'.format(label))
                fgnu.writelines('   set label "conv=".conv at 1.3,0.6\n')
                fgnu.writelines('   set label "merit=".merit at 1.3,0.4\n')
            fgnu.writelines('   p f1 u 1:2 w p pt 6 title "Iter ".n, TGT u 1:2 w l lc rgb "black" title "Target"\n'.format(label))
            if type_int == "bond" or type_int == "bend" or type_int == "dihedral":
                fgnu.writelines('   #unset label\n')
                fgnu.writelines('   #unset label\n')
            else:
                fgnu.writelines('   unset label\n')
                fgnu.writelines('   unset label\n')
            fgnu.writelines("}\n")

        # Run gnuplot to create the png file
        os.system("gnuplot {}".format(fnamegnu))

# =============================================================================
def calculate_merit_function(target_distributions, path_steps):

    for itarget, ivalue in target_distributions.items():
        type_int = ivalue["type"]
        if type_int != "non-bonded":
            continue

        dirs = sorted(glob.glob(os.path.join(path_steps, "step_*")))
        fmerit = open(os.path.join(os.getcwd(), "merit_" + itarget+".dat"), 'w')
        line_sum = ""
        for idir in dirs:
            fnamedist = os.path.join(idir, itarget+".dist.new")
            xdist = []
            ydist = []
            try:
                with open(fnamedist, 'r') as fdist:
                    lines = fdist.readlines()
                    for iline in lines:
                        if re.match("^#", iline):
                            continue
                        x, y, _ = iline.split()
                        xdist.append(x)
                        ydist.append(y)
                    a = np.array(ivalue["ydata"], dtype=np.float32)
                    b = np.array(ydist, dtype=np.float32)
                    aa = (a - b)*(a - b)
                    bb = b*b
                    sum1 = np.sum(aa)
                    sum2 = np.sum(bb)
                    line_sum += "{0:.2e}\n".format(sum1/sum2)
            except FileNotFoundError:
                continue
        fmerit.writelines(line_sum)
        fmerit.close()

    p = os.path.join(path_steps, "step_*")
    l = sorted(glob.glob(p))
    nsteps = len(l) - 1

    return nsteps

# =============================================================================
def prepare_individualpot_gnuplots(args, target_distributions, path_steps, terminal="png"):

    """
    Prepare gnu and png files for all distributions each step
    The png figures are located in the "path_steps+png" directory.

    Args:
        target_distributions: Dictionary with the info read from the VOTCA simulation
        path_steps: path where the steps_??? directories are located
        terminal: terminal line for gnuplot

    """

    if terminal == "png":
        terminal = "set terminal pngcairo size 550,400 enhanced font 'Verdana,10'\n" \
                   "set encoding iso_8859_1"
    else:
        terminal = terminal

    for itarget, ivalue in target_distributions.items():
        filename = ivalue['potfilename']
        type_int = ivalue["type"]
        label = os.path.basename(filename).split(".")[0]

        fnamegnu = "ibi_pot_"+label+".gnu"
        pattern_step = os.path.join(path_steps, "step_*")

        with open(fnamegnu, 'w') as fgnu:
            fgnu.writelines("reset\n")
            fgnu.writelines(terminal)
            fgnu.writelines("\n")
            fgnu.writelines("#Initial potential\n")
            fgnu.writelines('TGT="{}"\n'.format(filename))
            fgnu.writelines("\n")
            fgnu.writelines('DIR=system("ls -d {0:s} |egrep -v step_000")\n'.format(pattern_step))
            fgnu.writelines("\n")
            fgnu.writelines('set title "{} beads"\n'.format(label))
            if type_int == "bond":
                fgnu.writelines('set xlabel "r (nm)"\n')
                fgnu.writelines('set ylabel "E_{bond}/kbT"\n')
                fgnu.writelines('set format x "%.1f"\n')
                fgnu.writelines('set format y "%.1f"\n')
                fgnu.writelines('#set xrange [:]\n')
                fgnu.writelines('#set yrange [:]\n')
                fgnu.writelines('set grid\n')
            elif type_int == "bend":
                fgnu.writelines('set xlabel "{/Symbol q} (rad)"\n')
                fgnu.writelines('set ylabel "E_{angle}/kbT"\n')
                fgnu.writelines('set format x "%.1f"\n')
                fgnu.writelines('set format y "%.1f"\n')
                fgnu.writelines('#set xrange [:]\n')
                fgnu.writelines('#set yrange [:]\n')
                fgnu.writelines('set grid\n')
            elif type_int == "dihedral":
                fgnu.writelines('set xlabel "{/Symbol f} (rad)"\n')
                fgnu.writelines('set ylabel "E_{dih}/kbT"\n')
                fgnu.writelines('set format x "%.3f"\n')
                fgnu.writelines('set format y "%.3f"\n')
                fgnu.writelines('#set xrange [:]\n')
                fgnu.writelines('#set yrange [:]\n')
                fgnu.writelines('set grid\n')
            elif type_int == "non-bonded":
                fgnu.writelines('set xlabel "r (nm)"\n')
                fgnu.writelines('set ylabel "E_{{{0:s}}}(r)"\n'.format(label))
                fgnu.writelines('set format x "%.1f"\n')
                fgnu.writelines('set format y "%.1f"\n')
                fgnu.writelines('set xrange [0.4:1.2]\n')
                fgnu.writelines('set yrange [:2]\n')
                fgnu.writelines('set grid\n')
            fgnu.writelines('n=0\n')
            fgnu.writelines('do for [data in DIR] {\n')
            fgnu.writelines('   n=n+1\n')
            fgnu.writelines('   f1 = data."/{0:s}.pot.new"\n'.format(label))
            if type_int == "bond":
                fgnu.writelines('   set output sprintf("./png_pot/bond_{}_ibi%03d.png",n)\n'.format(label))
            if type_int == "bend":
                fgnu.writelines('   set output sprintf("./png_pot/bend_{}_ibi%03d.png",n)\n'.format(label))
            if type_int == "dihedral":
                fgnu.writelines('   set output sprintf("./png_pot/dihedral_{}_ibi%03d.png",n)\n'.format(label))
            if type_int == "non-bonded":
                fgnu.writelines('   set output sprintf("./png_pot/NB_{}_ibi%03d.png",n)\n'.format(label))
            fgnu.writelines('   p f1 u 1:($2/{0}) w p pt 6 title "Iter ".n, '
                            'TGT u 1:($2/{0}) w l lc rgb "black" title "Initial"\n'.format(args.kbT))
            fgnu.writelines("}\n")

        # Run gnuplot to create the png file
        os.system("gnuplot {}".format(fnamegnu))

# =============================================================================
def concatenate_images(nsteps, tdist, wkpath, rows=2, columns=3):

    from PIL import Image as PilImage

    # Check maximun number of figures
    if len(tdist) > 9:
        m = "\t\t More than 9 distributions are hard to visualize in one image.\n"
        m+= "\t\t The figures are not joined.\n"
        return m

    # Get all images from the directory and join the images
    join_filenames = []
    for istep in range(1, nsteps+1):
        images = []
        for ikey, ivalue in tdist.items():
            if ivalue["type"] != "non-bonded":
                fname = os.path.join(wkpath, "{0:s}_{1:s}_ibi{2:03d}.png".format(ivalue["type"], ikey, istep))
            else:
                fname = os.path.join(wkpath, "{0:s}_{1:s}_ibi{2:03d}.png".format("NB", ikey, istep))
            images.append(PilImage.open(fname))

        # Calculate size of the new image
        width = images[0].width * columns
        height = images[0].height * rows
        new_image = PilImage.new('RGB', (width, height))

        w = 0
        h = 0
        for irow in range(0, rows):
            for icol in range(0, columns):
                image_number = icol + (irow*columns)
                new_image.paste(images[image_number], (w, h))
                # Update width
                w += images[0].width
            # Update heigth
            w = 0
            h += images[0].height

        new_image.save("tmp_join_{0:03d}.png".format(istep))
        join_filenames.append("tmp_join_{0:03d}.png".format(istep))

    # Create the animated gif
    new_image = PilImage.new('RGB', (width, height))
    joinimages = []
    for i in join_filenames:
        joinimages.append(PilImage.open(i))
    fjoinname = "join_{}.gif".format(os.path.basename(wkpath))
    joinimages[0].save(fjoinname,
                       save_all=True, append_images=joinimages[1:], optimize=False, duration=400, loop=1)

    for itmp in glob.glob("tmp_join_*.png"):
        shutil.move(itmp, os.path.join(wkpath, itmp))

    m = "\t\t\t Joined figures are in the current directory."
    return m

# =============================================================================
def main_app():

    import polyanagro as pag

    # Parse arguments
    args = parse_arguments()
    # Setup log
    log = utils.init_logger("Output", fileoutput=args.log, append=False, inscreen=True)
    # Write header and arguments
    print_header(pag.version.__version__, log)

    # Prepare directory to store images
    d = os.path.join(os.getcwd(),"png_dist")
    shutil.rmtree(d, ignore_errors=True)
    os.mkdir("png_dist")

    # Distributions
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    d = os.path.join(args.path, "step_000")
    m = "\t\tRead target distributions from {} ({})".format(d, now)
    print(m) if log is None else log.info(m)
    target_distributions, type_distributions_dict = read_distribution_files(args)
    if target_distributions is None:
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        d = os.path.join(args.path, "step_000")
        m = "\t\tTarget distributions *.dist.tgt cannot be read from {} ({})".format(d, now)
        print(m) if log is None else log.error(m)
        exit()
    else:
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\tThe following distributions have been found:\n"
        m += "\t\t  Bond       : {}\n".format(len(type_distributions_dict["bond"]))
        for j in type_distributions_dict["bond"]:
            m += "\t\t      * {}\n".format(j)
        m += "\t\t  Bend       : {}\n".format(len(type_distributions_dict["bend"]))
        for j in type_distributions_dict["bend"]:
            m += "\t\t      * {}\n".format(j)
        m += "\t\t  Dihedral   : {}\n".format(len(type_distributions_dict["dihedral"]))
        for j in type_distributions_dict["dihedral"]:
            m += "\t\t      * {}\n".format(j)
        m += "\t\t  Non-bonded : {}\n".format(len(type_distributions_dict["non-bonded"]))
        for j in type_distributions_dict["non-bonded"]:
            m += "\t\t      * {}\n".format(j)
        print(m) if log is None else log.error(m)
        total_distributions = len(target_distributions)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tCalculating merit function for non-bonded iteractions ({})".format(now)
    print(m) if log is None else log.info(m)
    nsteps = calculate_merit_function(target_distributions, args.path)
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tPreparing gnuplot scripts and individual png graphs for distributions ({})".format(now)
    print(m) if log is None else log.info(m)
    prepare_individualdist_gnuplots(target_distributions, args.path)

    # Potentials
    dpot = os.path.join(os.getcwd(),"png_pot")
    shutil.rmtree(dpot, ignore_errors=True)
    os.mkdir("png_pot")

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tPreparing gnuplot scripts and individual png graphs for potentials ({})".format(now)
    print(m) if log is None else log.info(m)
    prepare_individualpot_gnuplots(args, target_distributions, args.path)

    if not args.tmpdir:
        shutil.rmtree(d, ignore_errors=True)
        shutil.rmtree(dpot, ignore_errors=True)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tCalculating merit function for non-bonded iteractions ({})".format(now)
    print(m) if log is None else log.info(m)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tJoining distributions to visualize ({})".format(now)
    print(m) if log is None else log.info(m)
    m = concatenate_images(nsteps, target_distributions, "./png_dist")
    m += concatenate_images(nsteps, target_distributions, "./png_pot")
    print(m) if log is None else log.info(m)

    m  = "\n\t\tRESULTS  : The png_dist and png_pot dirs contain the individual png graphs.\n"
    m2  = "\t\t           Compare new and target distributions.\n"
    m1 = "\t\t"+"-"*len(m)+"\n"
    print(m1+m+m2+m1) if log is None else log.info(m1+m+m2+m1)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tJob  Done at {} ============\n".format(now)
    print(m) if log is None else log.info(m)


# =============================================================================
if __name__ == "__main__":

    main_app()