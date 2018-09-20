import shutil
import subprocess
import os

### make directories ###
def code_mkdirs():
    f_dirname = open("code_output/dirnames")
    dirnames = f_dirname.readlines()
    f_dirname.close()
    for i in range(0,len(dirnames)):
        dirnames[i] = dirnames[i].rstrip()
        if not os.path.isdir(dirnames[i]):
            os.mkdir(dirnames[i])
            print("make dir "+dirnames[i])
        else:
            print("Already exist dir "+dirnames[i])
    print("FINISH mkdirs")

### binfile ###
def code_binfile():
    workdir = "binfile"
    shutil.copy("code_output/merge_4ES.f90", "binfile/")
    subprocess.check_call(["ifort", "merge_4ES.f90"],
                          cwd=workdir)
    subprocess.check_call("./a.out",
                          cwd=workdir)
    print("FINISH code_binfile")
### 2d/xy ###
def code_2dxy():
    workdir = "2d/xy"
    shutil.copy("code_output/draw_xy.py","2d/xy/")
    subprocess.check_call(["python", "draw_xy.py"],
                          cwd=workdir)
    print("FINISH code_2dxy")
### 2d/yz ###
def code_2dyz():
    workdir = "2d/yz"
    shutil.copy("code_output/draw_yz.py","2d/yz/")
    subprocess.check_call(["python", "draw_yz.py"],
                          cwd=workdir)
    print("FINISH code_2dyz")
### 3dbin ###
def code_3dbin():
    workdir = "3dbin"
    shutil.copy("code_output/convert.f90", "3dbin/")
    subprocess.check_call(["ifort", "convert.f90"],
                          cwd=workdir)
    subprocess.check_call("./a.out",
                          cwd=workdir)
    print("FINISH code_3dbin")

### r3dbline_calculation ###
def code_r3dbline_calc():
    workdir = "code_output/code_bline/"
    subprocess.check_call(["make", "clean"],
                          cwd=workdir)
    subprocess.check_call("make",
                          cwd=workdir)
    shutil.copy("code_output/code_bline/bline_xy", "r3dbline/")
    shutil.copy("code_output/code_bline/NAMELIST", "r3dbline/")
    shutil.copy("code_output/code_bline/Param_bline", "r3dbline/")
    workdir = "r3dbline/"
    subprocess.check_call("./bline_xy",
                          cwd=workdir)
    print("FINISH code_r3dbline_calc")

### r3dbline_plotmaps ###
def code_r3dbline_plot():
    ### copy script files
    shutil.copy("code_output/code_pltbline/plotbyeq0.py",
                "r3dbline/file_output/picture/byeq0/")
    shutil.copy("code_output/code_pltbline/draw_connect.py",
                "r3dbline/file_output/picture/connect/")
    shutil.copy("code_output/code_pltbline/mkfile.f90",
                "r3dbline/file_output/picture/foot/")
    shutil.copy("code_output/code_pltbline/drawfoot.py",
                "r3dbline/file_output/picture/foot/")
    shutil.copy("code_output/code_pltbline/detectHTR.f90",
                "r3dbline/file_output/picture/hightwist/")
    shutil.copy("code_output/code_pltbline/drawjjyz.py",
                "r3dbline/file_output/picture/jjyz/")
    shutil.copy("code_output/code_pltbline/calc_fairec.f90",
                "r3dbline/file_output/picture/recflux/")
    shutil.copy("code_output/code_pltbline/output_recflux.py",
                "r3dbline/file_output/picture/recflux/")
    shutil.copy("code_output/code_pltbline/draw_ts.py",
                "r3dbline/file_output/picture/timeslice_FR/")
    shutil.copy("code_output/code_pltbline/drawtw.py",
                "r3dbline/file_output/picture/tw/")
    ### plot map of by.eq.0 ###
    workdir = "r3dbline/file_output/picture/byeq0/"
    subprocess.check_call(["python","plotbyeq0.py"],
                          cwd=workdir)
    ### calculate footpoints info & plot footpoints map ###
    workdir = "r3dbline/file_output/picture/foot/"
    subprocess.check_call(["ifort","mkfile.f90"],
                          cwd=workdir)
    subprocess.check_call("./a.out",
                          cwd=workdir)
    subprocess.check_call(["python","drawfoot.py"],
                          cwd=workdir)
    ### plot connection map ###
    # workdir = "r3dbline/file_output/picture/connect/"
    # subprocess.check_call(["python","draw_connect.py"],
    #                       cwd=workdir)
    ### integrate high twist region ###
    # workdir = "r3dbline/file_output/picture/hightwist/"
    # subprocess.check_call(["ifort","detectHTR.f90"],
    #                       cwd=workdir)
    # subprocess.check_call("./a.out",
    #                       cwd=workdir)
    ### draw jj on yz plane ###
    workdir = "r3dbline/file_output/picture/jjyz/"
    subprocess.check_call(["python","drawjjyz.py"],
                          cwd=workdir)
    ### integrate reconnected flux and plot ###
    workdir = "r3dbline/file_output/picture/recflux/"
    subprocess.check_call(["ifort","calc_fairec.f90"],
                          cwd=workdir)
    subprocess.check_call("./a.out",
                          cwd=workdir)
    subprocess.check_call(["python","output_recflux.py"],
                          cwd=workdir)
    ### plot time slice on yz plane ###
    workdir = "r3dbline/file_output/picture/timeslice_FR/"
    subprocess.check_call(["python","draw_ts.py"],
                          cwd=workdir)
    ### make twist map ###
    workdir = "r3dbline/file_output/picture/tw/"
    subprocess.check_call(["python","drawtw.py"],
                          cwd=workdir)

### calculate and output kappa ###
def code_calckappa():
    ### copy script files
    shutil.copy("code_output/kappa.py",
                "calc_kappa/")
    ### calculate and output kappa ###
    workdir = "calc_kappa/"
    subprocess.check_call(["python","kappa.py"],
                          cwd=workdir)

### make vapor .vdf files ###
def code_vapor_makefile():
    ### copy script files
    shutil.copy("code_output/vdf_cre.sh",
                "vapor/")
    ### execute script file for make .vdf files ###
    workdir = "vapor/"
    subprocess.check_call("./vdf_cre.sh",
                          cwd=workdir)

### calculate and output kappa2 ###
def code_kappa2():
    ### copy script files
    shutil.copy("code_output/kappa2.py",
                "kappa2/")
    ### calculate and output kappa ###
    workdir = "kappa2/"
    subprocess.check_call(["python","kappa2.py"],
                          cwd=workdir)
