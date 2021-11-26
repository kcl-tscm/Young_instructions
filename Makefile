#########################################################################
# Makefile for cpmd.x (plane wave electronic calculation)
# Configuration: LINUX-X86_64-INTEL-MPI
# Creation of Makefile: Nov 26 2021
# on Linux login01 3.10.0-1127.13.1.el7.x86_64 #1 SMP Fri Jun 12 14:34:17 EDT 2020 x86_64 x86_64 x86_64 GNU/Linux
# Author: mmm0666
#
# To compile: go into the DEST directory (if using DEST) or stay here 
# and type make.
#########################################################################
#
SHELL = /bin/sh

#########################################################################
# project-independent stuff

CPMDROOT = /home/mmm0666/CPMD_4.1_SVN_QMMM/CPMD
DEST     = /home/mmm0666/CPMD_4.1_SVN_QMMM/CPMD/
SRCDIR   = ${CPMDROOT}/src
MODDIR   = ${CPMDROOT}/modules

BINDIR   = ${DEST}/bin
OBJDIR   = ${DEST}/obj
LIBDIR   = ${DEST}/lib

MAKEFILE = ${DEST}/Makefile

TARGET        = $(BINDIR)/cpmd.x
CPMD_LIB      = $(LIBDIR)/libcpmd.a
GROMOS_LIB    = $(LIBDIR)/libgromos.a
INTERFACE_LIB = $(LIBDIR)/libinterface.a

.SUFFIXES: .F90 .f90 .c .o

#########################################################################
# platform-dependent stuff
# this section is built by the configure.sh script: no manual editing
# should be required.
#########################################################################
FFLAGS = -O2 -pc64 -funroll-loops -qopenmp -I${SRCDIR} -I${OBJDIR}
##standard MKL parallel
#LFLAGS = -static-intel -mkl=parallel
##standard MKL sequential
LFLAGS = -static-intel -mkl=parallel
##Personal BLAS/LAPACK
#LFLAGS = -L/home/boero/lib/LAPACK -llapack-amd -lblas-amd
CFLAGS = -O2 -Wall -m64 -I${SRCDIR}
NVCCFLAGS =  -I${SRCDIR}
CPP = fpp -P -C
CPPFLAGS = -D__Linux -D__INTEL -DLINUX_IFC -D__HAS_FFT_DEFAULT -D__HASNT_F08_ISO_FORTRAN_ENV \
           -D__HPC -D__HAS_SIZEOF -D__PARALLEL  -D__GROMOS  \
           -I${SRCDIR} -D'SVN_REV="r$(shell svnversion -n /home/mmm0666/CPMD_4.1_SVN_QMMM/CPMD)"'
NOOPT_FLAG = 
CPPFLAGS_GROMOS = -DEWALD -DEWATCUT -DHAT_SHAPE -DUNPACKED_GRID  
#FFLAGS_GROMOS = -fixed $(FFLAGS) -I${MODDIR}
#FFLAGS_GROMOS_MODULES = $(FFLAGS) -I${MODDIR}
FFLAGS_GROMOS = -fixed -O2 -pc64 -funroll-loops -I${SRCDIR} -I${OBJDIR} -I${MODDIR}
FFLAGS_GROMOS_MODULES = -O2 -pc64 -funroll-loops -I${SRCDIR} -I${OBJDIR} -I${MODDIR}
FFLAGS += -I/home/mmm0666/CPMD_4.1_SVN_QMMM/CPMD/modules
CC = mpicc
FC = mpif90
LD = mpif90
NVCC = 
AR = xiar -r
#AR = ar -r
RANLIB = ranlib
#########################################################################
# Personal Configuration
# My_Conf: 
# All arguments:  LINUX-X86_64-INTEL-MPI -qmmm
CONFARGS = LINUX-X86_64-INTEL-MPI -qmmm
#########################################################################
#########################################################################
# End of Personal Configuration
#########################################################################
CFGDEST = /home/mmm0666/CPMD_4.1_SVN_QMMM/CPMD//obj
CFGMACH = LINUX-X86_64-INTEL-MPI
CFGQMMM = -qmmm
################################################################################
# files 
# load SRC_ variables from SOURCES file: correct module compilation order
include $(SRCDIR)/SOURCES

OBJ_MODULES = $(SRC_MODULES:%.F90=%.o)
OBJ_F90     = $(SRC_F90:%.F90=%.o)
OBJ_CC      = $(SRC_CC:%.c=%.o)
ifneq ($(NVCC),)
OBJ_CU      = $(SRC_CU:%.cu=%.o)
else
OBJ_CU      = 
endif
OBJ_CPMD     = $(OBJ_MODULES) $(OBJ_F90) $(OBJ_CC) $(OBJ_CU)

# load QMMM SRC_ variables from SOURCES file
include $(MODDIR)/QMMM_SOURCES

OBJECTS_MODULES_GROMOS = $(GROMOS_MODULES:%.F90=%.o) 
OBJECTS_GROMOS    = $(GROMOS_SRC:%.F=%.o)
OBJECTS_INTERFACE = $(INTERFACE_SRC:%.F=%.o)

################################################################################
# Explicit Rules
################################################################################


default, all:
	+$(MAKE) -C $(OBJDIR) -f $(MAKEFILE) $(TARGET)

lib:
	+$(MAKE) -C $(OBJDIR) -f $(MAKEFILE) $(CPMD_LIB)
	+$(MAKE) -C $(OBJDIR) -f $(MAKEFILE) $(GROMOS_LIB)
	+$(MAKE) -C $(OBJDIR) -f $(MAKEFILE) $(INTERFACE_LIB)

distclean:  realclean
	rm -f  $(TARGET)
	rm -rf $(DEST)/Makefile

realclean:  clean
	rm -rf $(OBJDIR) $(BINDIR) $(LIBDIR)

clean:  purge
	rm -f $(CPMD_LIB) $(GROMOS_LIB) $(INTERFACE_LIB) $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/*.f90 $(OBJDIR)/*.c $(OBJDIR)/*.cu
	rm -f $(OBJDIR)/Gromos/*.o $(OBJDIR)/Gromos/*.f*
	rm -f $(OBJDIR)/MM_Interface/*.o $(OBJDIR)/MM_Interface/*.f
	rm -f $(OBJDIR)/IPhigenie_Interface/*.{o,f}

purge: 
	rm -f $(SRCDIR)/*~

reconfigure:
	$(shell cd $(CPMDROOT); ./configure.sh $(CONFARGS) )

%.o:    $(SRCDIR)/%.c
	( cd $(SRCDIR); cp $(<F) $(OBJDIR)/$(<F); cd $(OBJDIR) )
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $(<F)

ifneq ($(NVCC),)
%.o:    $(SRCDIR)/%.cu
	( cd $(SRCDIR); cp $(<F) $(OBJDIR)/$(<F); cd $(OBJDIR) )
	$(NVCC) -c $(NVCCFLAGS) $(CPPFLAGS) -o $@ $(<F)
endif

ifneq ($(CPP),)
%.f90:  $(SRCDIR)/%.F90
	( cd $(SRCDIR); $(CPP) $(CPPFLAGS) $(<F) $(OBJDIR)/$(@F); cd $(OBJDIR) )

%.mod.o: %.mod.f90
	$(FC) -c $(FFLAGS) -o $@ $(@F:.o=.f90)

%.o:    $(OBJ_MODULES) %.f90 
	$(FC) -c $(FFLAGS) -o $@ $(@F:.o=.f90)
else
%.mod.o: $(SRCDIR)/%.mod.F90
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $@ $(SRCDIR)/$(@F:.o=.F90)

%.o:    $(OBJ_MODULES) %.f90 
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $@ $(SRCDIR)/$(@F:.o=.F90)
endif
$(TARGET): $(CPMD_LIB) $(GROMOS_LIB) $(INTERFACE_LIB) timetag.o cpmd.o
	$(LD) $(FFLAGS) -o $(TARGET) timetag.o cpmd.o $(CPMD_LIB) $(GROMOS_LIB) $(INTERFACE_LIB) $(LFLAGS)
	@ls -l $(TARGET)
	@echo "Compilation done."

$(CPMD_LIB): $(OBJ_CPMD)
	$(AR) $(CPMD_LIB) $(OBJ_CPMD)
	$(RANLIB) $(CPMD_LIB)

################################################################################
#       SPECIAL RULES ONLY
################################################################################
# IRAT module is compiled separately to have -D_IRAT_
irat.mod.o: $(SRCDIR)/irat.mod.F90
ifneq ($(CPP),) 
	$(CPP) $(CPPFLAGS) -D_IRAT_=2 $< $(@F:.o=.f90)
	$(FC) -c $(FFLAGS) -o $@ $(@F:.o=.f90)
else
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -D_IRAT_=2 -o $@ $(SRCDIR)/$(@F:.o=.F90)
endif

# timetag is compiled separately to have the compilation date always up2date
.PHONY: timetag.o
timetag.o:
ifneq ($(CPP),) 
	$(CPP) $(CPPFLAGS) $(SRCDIR)/timetag.F90 timetag.f90
	$(FC) -c $(FFLAGS) -o timetag.o timetag.f90
else
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o timetag.o $(SRCDIR)/timetag.F90 
endif

.SECONDARY: $(SRC_ALL:%.F90=%.f90) cpmd.f90 # do not delete intermediate .f90 files
.SECONDARY: $(GROMOS_SRC:%.F=%.f) $(INTERFACE_SRC:%.F=%.f)
.SECONDARY: $(GROMOS_MODULES:%.mod.F90=%.mod.f90) 

################################################################################
#       QM/MM rules
################################################################################
# TODO write special rules for forcematch_ only
$(OBJ_MODULES): $(OBJECTS_MODULES_GROMOS)

ifneq ($(CPP),)
$(OBJECTS_GROMOS:.o=.f): 
	rm -f $@
	$(CPP) $(CPPFLAGS_GROMOS) $(MODDIR)/$(@:.f=.F) $@

$(OBJECTS_MODULES_GROMOS:.o=.f90):
	rm -f $@
	$(CPP) $(CPPFLAGS_GROMOS) $(MODDIR)/$(@:.f90=.F90) $@

$(OBJECTS_MODULES_GROMOS): $(OBJECTS_MODULES_GROMOS:.o=.f90)

	$(FC) -c $(FFLAGS_GROMOS_MODULES) -I$(MODDIR)/Gromos/ -I$(MODDIR)/MM_Interface/ -o $@ $(@:.o=.f90) 

$(OBJECTS_GROMOS): $(OBJECTS_MODULES_GROMOS)
	$(FC) -c $(FFLAGS_GROMOS) -I$(MODDIR)/Gromos/ -I$(MODDIR)/MM_Interface/ -o $@ $(@:.o=.f)

else
$(OBJECTS_MODULES_GROMOS):
	$(FC) -c $(FFLAGS_GROMOS_MODULES) $(CPPFLAGS_GROMOS) -I$(MODDIR)/Gromos/ -I$(MODDIR)/MM_Interface/ -o $@ $(MODDIR)/$(@:.o=.F90) 

$(OBJECTS_GROMOS): $(OBJECTS_MODULES_GROMOS) 
	$(FC) -c $(FFLAGS_GROMOS) $(CPPFLAGS_GROMOS) -I$(MODDIR)/Gromos/ -I$(MODDIR)/MM_Interface/ -o $@ $(MODDIR)/$(@:.o=.F) 
endif

$(GROMOS_LIB): $(OBJECTS_MODULES_GROMOS) $(OBJECTS_GROMOS)
	$(AR) $(GROMOS_LIB) $(OBJECTS_MODULES_GROMOS) $(OBJECTS_GROMOS)
	$(RANLIB) $(GROMOS_LIB)

ifneq ($(CPP),)
$(OBJECTS_INTERFACE:.o=.f): 
	rm -f $@
	$(CPP) $(CPPFLAGS) $(MODDIR)/$(@:.f=.F) $@

$(OBJECTS_INTERFACE): 
	$(FC) -c $(FFLAGS_GROMOS) -I$(MODDIR)/Gromos/ -I$(MODDIR)/MM_Interface/ $< -o $@
else
$(OBJECTS_INTERFACE): 
	$(FC) -c $(FFLAGS_GROMOS) $(CPPFLAGS_GROMOS) -I$(MODDIR)/Gromos/ -I$(MODDIR)/MM_Interface/ -o $@ $(MODDIR)/$(@:.o=.F)
endif

$(INTERFACE_LIB): $(OBJECTS_INTERFACE)
	$(AR) $(INTERFACE_LIB) $(OBJECTS_INTERFACE)
	$(RANLIB) $(INTERFACE_LIB)
################################################################################
# Module dependencies
# 
aainit_utils.mod:aainit_utils.mod.o
	@true
aavan.mod:      aavan.mod.o
	@true
adapttol_utils.mod:adapttol_utils.mod.o
	@true
adat.mod:       adat.mod.o
	@true
adjmu_utils.mod:adjmu_utils.mod.o
	@true
afbdr_utils.mod:afbdr_utils.mod.o
	@true
ainitwf_utils.mod:ainitwf_utils.mod.o
	@true
anderson_utils.mod:anderson_utils.mod.o
	@true
andp.mod:       andp.mod.o
	@true
andr.mod:       andr.mod.o
	@true
anneal_utils.mod:anneal_utils.mod.o
	@true
array_utils.mod:array_utils.mod.o
	@true
atimesmod.mod:  atimes.mod.o
	@true
atimes_utils.mod:atimes_utils.mod.o
	@true
atomc_utils.mod:atomc_utils.mod.o
	@true
atom.mod:       atom.mod.o
	@true
atoms_utils.mod:atoms_utils.mod.o
	@true
atomwf_utils.mod:atomwf_utils.mod.o
	@true
atrho_utils.mod:atrho_utils.mod.o
	@true
atwf.mod:       atwf.mod.o
	@true
augchg_utils.mod:augchg_utils.mod.o
	@true
azzero_utils.mod:azzero_utils.mod.o
	@true
bc.mod:         bc.mod.o
	@true
benc.mod:       benc.mod.o
	@true
bessm_utils.mod:bessm_utils.mod.o
	@true
bicanonicalCalculationConfig.mod:bicanonicalCalculationConfig.mod.o
	@true
bicanonicalCalculation.mod:bicanonicalCalculation.mod.o
	@true
bicanonicalConfig.mod:bicanonicalConfig.mod.o
	@true
bicanonicalCpmd.mod:bicanonicalCpmd.mod.o
	@true
bicanonicalInputReader.mod:bicanonicalInputReader.mod.o
	@true
bogol_utils.mod:bogol_utils.mod.o
	@true
box_boundary_utils.mod:box_boundary_utils.mod.o
	@true
broyden_utils.mod:broyden_utils.mod.o
	@true
broy.mod:       broy.mod.o
	@true
bs_forces_diag_utils.mod:bs_forces_diag_utils.mod.o
	@true
bswfo_utils.mod:bswfo_utils.mod.o
	@true
bsym.mod:       bsym.mod.o
	@true
bsympnt.mod:    bsympnt.mod.o
	@true
calc_alm_utils.mod:calc_alm_utils.mod.o
	@true
calc_pij_utils.mod:calc_pij_utils.mod.o
	@true
canon_utils.mod:canon_utils.mod.o
	@true
cdftmod.mod:    cdft.mod.o
	@true
cdft_utils.mod: cdft_utils.mod.o
	@true
cell.mod:       cell.mod.o
	@true
chain_dr_utils.mod:chain_dr_utils.mod.o
	@true
chksym_utils.mod:chksym_utils.mod.o
	@true
clas_force_utils.mod:clas_force_utils.mod.o
	@true
clas.mod:       clas.mod.o
	@true
clinbcg_utils.mod:clinbcg_utils.mod.o
	@true
cl_init_utils.mod:cl_init_utils.mod.o
	@true
cmaos_utils.mod:cmaos_utils.mod.o
	@true
cnst_dyn.mod:   cnst_dyn.mod.o
	@true
cnstfc_utils.mod:cnstfc_utils.mod.o
	@true
cnst.mod:       cnst.mod.o
	@true
cnstpr_utils.mod:cnstpr_utils.mod.o
	@true
cofor_utils.mod:cofor_utils.mod.o
	@true
comvelmod.mod:  comvel.mod.o
	@true
comvel_utils.mod:comvel_utils.mod.o
	@true
conduct_utils.mod:conduct_utils.mod.o
	@true
condu.mod:      condu.mod.o
	@true
coninp_utils.mod:coninp_utils.mod.o
	@true
constr_utils.mod:constr_utils.mod.o
	@true
control_bcast_utils.mod:control_bcast_utils.mod.o
	@true
control_def_utils.mod:control_def_utils.mod.o
	@true
control_pri_utils.mod:control_pri_utils.mod.o
	@true
control_test_utils.mod:control_test_utils.mod.o
	@true
control_utils.mod:control_utils.mod.o
	@true
conv.mod:       conv.mod.o
	@true
coor.mod:       coor.mod.o
	@true
copot_utils.mod:copot_utils.mod.o
	@true
corec_utils.mod:corec_utils.mod.o
	@true
cores.mod:      cores.mod.o
	@true
core_spect_utils.mod:core_spect_utils.mod.o
	@true
cotr.mod:       cotr.mod.o
	@true
cp_cuda_types.mod:cp_cuda_types.mod.o
	@true
cp_cuda_utils.mod:cp_cuda_utils.mod.o
	@true
cp_cudensity_utils.mod:cp_cudensity_utils.mod.o
	@true
cp_cufft_types.mod:cp_cufft_types.mod.o
	@true
cp_cufft_utils.mod:cp_cufft_utils.mod.o
	@true
cp_cuortho_types.mod:cp_cuortho_types.mod.o
	@true
cp_cuortho_utils.mod:cp_cuortho_utils.mod.o
	@true
cp_curho_types.mod:cp_curho_types.mod.o
	@true
cp_curho_utils.mod:cp_curho_utils.mod.o
	@true
cp_cuvpsi_types.mod:cp_cuvpsi_types.mod.o
	@true
cp_cuvpsi_utils.mod:cp_cuvpsi_utils.mod.o
	@true
cp_cuwfn_types.mod:cp_cuwfn_types.mod.o
	@true
cp_cuwfn_utils.mod:cp_cuwfn_utils.mod.o
	@true
cp_dgga_correlation_utils.mod:cp_dgga_correlation_utils.mod.o
	@true
cp_dgga_exchange_utils.mod:cp_dgga_exchange_utils.mod.o
	@true
cp_dxc_driver.mod:cp_dxc_driver.mod.o
	@true
cpfunc_types.mod:cpfunc_types.mod.o
	@true
cp_gga_correlation_utils.mod:cp_gga_correlation_utils.mod.o
	@true
cp_gga_exchange_utils.mod:cp_gga_exchange_utils.mod.o
	@true
cp_grp_utils.mod:cp_grp_utils.mod.o
	@true
cp_ieee_interface.mod:cp_ieee_interface.mod.o
	@true
cp_lda_correlation_utils.mod:cp_lda_correlation_utils.mod.o
	@true
cp_lda_exchange_utils.mod:cp_lda_exchange_utils.mod.o
	@true
cplngsmod.mod:  cplngs.mod.o
	@true
cplngs_utils.mod:cplngs_utils.mod.o
	@true
cp_mgga_correlation_utils.mod:cp_mgga_correlation_utils.mod.o
	@true
cp_mgga_exchange_utils.mod:cp_mgga_exchange_utils.mod.o
	@true
cppt.mod:       cppt.mod.o
	@true
cp_xc_driver.mod:cp_xc_driver.mod.o
	@true
cp_xc_utils.mod:cp_xc_utils.mod.o
	@true
crotwf_utils.mod:crotwf_utils.mod.o
	@true
csize_utils.mod:csize_utils.mod.o
	@true
csmat_utils.mod:csmat_utils.mod.o
	@true
cublas_interfaces.mod:cublas_interfaces.mod.o
	@true
cublas_types.mod:cublas_types.mod.o
	@true
cublas_utils.mod:cublas_utils.mod.o
	@true
cuda_interfaces.mod:cuda_interfaces.mod.o
	@true
cuda_types.mod: cuda_types.mod.o
	@true
cuda_utils.mod: cuda_utils.mod.o
	@true
cufft_interfaces.mod:cufft_interfaces.mod.o
	@true
cufft_types.mod:cufft_types.mod.o
	@true
cufft_utils.mod:cufft_utils.mod.o
	@true
cusolver_interfaces.mod:cusolver_interfaces.mod.o
	@true
cusolver_types.mod:cusolver_types.mod.o
	@true
cusolver_utils.mod:cusolver_utils.mod.o
	@true
cuuser_interfaces.mod:cuuser_interfaces.mod.o
	@true
cuuser_utils.mod:cuuser_utils.mod.o
	@true
cvan.mod:       cvan.mod.o
	@true
davidson_utils.mod:davidson_utils.mod.o
	@true
dcacp_utils.mod:dcacp_utils.mod.o
	@true
dd_functionals_utils.mod:dd_functionals_utils.mod.o
	@true
ddip.mod:       ddip.mod.o
	@true
ddipo_utils.mod:ddipo_utils.mod.o
	@true
dd_xc_ana_utils.mod:dd_xc_ana_utils.mod.o
	@true
dd_xc_utils.mod:dd_xc_utils.mod.o
	@true
debfor_utils.mod:debfor_utils.mod.o
	@true
density_functionals_utils.mod:density_functionals_utils.mod.o
	@true
density_utils.mod:density_utils.mod.o
	@true
densrd_utils.mod:densrd_utils.mod.o
	@true
densto_utils.mod:densto_utils.mod.o
	@true
deort_utils.mod:deort_utils.mod.o
	@true
detdof_utils.mod:detdof_utils.mod.o
	@true
detsp_utils.mod:detsp_utils.mod.o
	@true
dftin_utils.mod:dftin_utils.mod.o
	@true
dginit_utils.mod:dginit_utils.mod.o
	@true
dg.mod:         dg.mod.o
	@true
difrho_utils.mod:difrho_utils.mod.o
	@true
dipomod.mod:    dipo.mod.o
	@true
dipo_utils.mod: dipo_utils.mod.o
	@true
disortho_utils.mod:disortho_utils.mod.o
	@true
dispp_utils.mod:dispp_utils.mod.o
	@true
dist_friesner_utils.mod:dist_friesner_utils.mod.o
	@true
dist_prowfn_utils.mod:dist_prowfn_utils.mod.o
	@true
d_mat_p_utils.mod:d_mat_p_utils.mod.o
	@true
dnlpdk_p_utils.mod:dnlpdk_p_utils.mod.o
	@true
domdr_utils.mod:domdr_utils.mod.o
	@true
do_perturbation_p_utils.mod:do_perturbation_p_utils.mod.o
	@true
dotp_utils.mod: dotp_utils.mod.o
	@true
dpot.mod:       dpot.mod.o
	@true
dqgalloc_utils.mod:dqgalloc_utils.mod.o
	@true
dqvan2_utils.mod:dqvan2_utils.mod.o
	@true
drhov_utils.mod:drhov_utils.mod.o
	@true
dum2_utils.mod: dum2_utils.mod.o
	@true
dylmr_utils.mod:dylmr_utils.mod.o
	@true
dynit_utils.mod:dynit_utils.mod.o
	@true
eam.mod:        eam.mod.o
	@true
eam_pot_utils.mod:eam_pot_utils.mod.o
	@true
eextern_utils.mod:eextern_utils.mod.o
	@true
efield_utils.mod:efield_utils.mod.o
	@true
efld.mod:       efld.mod.o
	@true
egointer_utils.mod:egointer_utils.mod.o
	@true
ehpsi_utils.mod:ehpsi_utils.mod.o
	@true
ehrenfest_utils.mod:ehrenfest_utils.mod.o
	@true
eicalc_utils.mod:eicalc_utils.mod.o
	@true
eigensystem_p_utils.mod:eigensystem_p_utils.mod.o
	@true
eind_ii_utils.mod:eind_ii_utils.mod.o
	@true
eind_loc_utils.mod:eind_loc_utils.mod.o
	@true
eind_nl_utils.mod:eind_nl_utils.mod.o
	@true
ekinpp_utils.mod:ekinpp_utils.mod.o
	@true
elct2.mod:      elct2.mod.o
	@true
elct.mod:       elct.mod.o
	@true
elf_utils.mod:  elf_utils.mod.o
	@true
elstpo_utils.mod:elstpo_utils.mod.o
	@true
empf.mod:       empf.mod.o
	@true
empfor_utils.mod:empfor_utils.mod.o
	@true
enbandpri_utils.mod:enbandpri_utils.mod.o
	@true
ener.mod:       ener.mod.o
	@true
enosmove_utils.mod:enosmove_utils.mod.o
	@true
envir_utils.mod:envir_utils.mod.o
	@true
envj.mod:       envj.mod.o
	@true
epot_types.mod: epot_types.mod.o
	@true
epr_current_p_utils.mod:epr_current_p_utils.mod.o
	@true
epr_dv0_utils.mod:epr_dv0_utils.mod.o
	@true
epr_efg_utils.mod:epr_efg_utils.mod.o
	@true
epr_hyp_utils.mod:epr_hyp_utils.mod.o
	@true
epr_p_utils.mod:epr_p_utils.mod.o
	@true
epr_util_p_utils.mod:epr_util_p_utils.mod.o
	@true
error_handling.mod:error_handling.mod.o
	@true
espchg_utils.mod:espchg_utils.mod.o
	@true
evirial_utils.mod:evirial_utils.mod.o
	@true
exdipo_utils.mod:exdipo_utils.mod.o
	@true
exterp_utils.mod:exterp_utils.mod.o
	@true
extpotmod.mod:  extpot.mod.o
	@true
extrap_utils.mod:extrap_utils.mod.o
	@true
fcas.mod:       fcas.mod.o
	@true
ffsum_utils.mod:ffsum_utils.mod.o
	@true
fftchk_utils.mod:fftchk_utils.mod.o
	@true
fftcu_methods.mod:fftcu_methods.mod.o
	@true
fftmain_utils.mod:fftmain_utils.mod.o
	@true
fft_maxfft.mod: fft_maxfft.mod.o
	@true
fft.mod:        fft.mod.o
	@true
fftnew_utils.mod:fftnew_utils.mod.o
	@true
fftprp_utils.mod:fftprp_utils.mod.o
	@true
fft_utils.mod:  fft_utils.mod.o
	@true
fftutil_utils.mod:fftutil_utils.mod.o
	@true
fharm_utils.mod:fharm_utils.mod.o
	@true
fileopenmod.mod:fileopen.mod.o
	@true
fileopen_utils.mod:fileopen_utils.mod.o
	@true
fillc_utils.mod:fillc_utils.mod.o
	@true
filnmod.mod:    filn.mod.o
	@true
finalp_utils.mod:finalp_utils.mod.o
	@true
fint.mod:       fint.mod.o
	@true
fitpack_utils.mod:fitpack_utils.mod.o
	@true
fixcom_utils.mod:fixcom_utils.mod.o
	@true
fm_cnst.mod:    fm_cnst.mod.o
	@true
fnlalloc_utils.mod:fnlalloc_utils.mod.o
	@true
fnonloc_p_utils.mod:fnonloc_p_utils.mod.o
	@true
fnonloc_utils.mod:fnonloc_utils.mod.o
	@true
forcedr_driver.mod:forcedr_driver.mod.o
	@true
forcedr_utils.mod:forcedr_utils.mod.o
	@true
forcematch_kfit_utils.mod:forcematch_kfit_utils.mod.o
	@true
forcematch.mod: forcematch.mod.o
	@true
forcematch_qfit_utils.mod:forcematch_qfit_utils.mod.o
	@true
forcematch_utils.mod:forcematch_utils.mod.o
	@true
forcep_utils.mod:forcep_utils.mod.o
	@true
forces_diag_utils.mod:forces_diag_utils.mod.o
	@true
forces_driver.mod:forces_driver.mod.o
	@true
forces_prop_utils.mod:forces_prop_utils.mod.o
	@true
forces_p_utils.mod:forces_p_utils.mod.o
	@true
forces_utils.mod:forces_utils.mod.o
	@true
formf_utils.mod:formf_utils.mod.o
	@true
freqs_utils.mod:freqs_utils.mod.o
	@true
friesner_c_p_utils.mod:friesner_c_p_utils.mod.o
	@true
friesner_c_utils.mod:friesner_c_utils.mod.o
	@true
friesner_utils.mod:friesner_utils.mod.o
	@true
frsblk_c_utils.mod:frsblk_c_utils.mod.o
	@true
frsblk_utils.mod:frsblk_utils.mod.o
	@true
fstart_utils.mod:fstart_utils.mod.o
	@true
fukui_p_utils.mod:fukui_p_utils.mod.o
	@true
func.mod:       func.mod.o
	@true
functionals_utils.mod:functionals_utils.mod.o
	@true
fusion_utils.mod:fusion_utils.mod.o
	@true
gcener_utils.mod:gcener_utils.mod.o
	@true
gcxctbl_utils.mod:gcxctbl_utils.mod.o
	@true
genxc_utils.mod:genxc_utils.mod.o
	@true
geofile_utils.mod:geofile_utils.mod.o
	@true
geq0mod.mod:    geq0.mod.o
	@true
getcor_utils.mod:getcor_utils.mod.o
	@true
getfnm_utils.mod:getfnm_utils.mod.o
	@true
getfu_utils.mod:getfu_utils.mod.o
	@true
getgyr_utils.mod:getgyr_utils.mod.o
	@true
gettrans_utils.mod:gettrans_utils.mod.o
	@true
gfft_utils.mod: gfft_utils.mod.o
	@true
ghermit_utils.mod:ghermit_utils.mod.o
	@true
glemod.mod:     gle.mod.o
	@true
gle_utils.mod:  gle_utils.mod.o
	@true
global_utils.mod:global_utils.mod.o
	@true
g_loc_dr_utils.mod:g_loc_dr_utils.mod.o
	@true
g_loc_exp_sum_utils.mod:g_loc_exp_sum_utils.mod.o
	@true
g_loc.mod:      g_loc.mod.o
	@true
g_loc_opeigr_utils.mod:g_loc_opeigr_utils.mod.o
	@true
g_loc_optim_utils.mod:g_loc_optim_utils.mod.o
	@true
g_loc_realspace_utils.mod:g_loc_realspace_utils.mod.o
	@true
g_loc_spread_ide_utils.mod:g_loc_spread_ide_utils.mod.o
	@true
g_loc_spread_sum_utils.mod:g_loc_spread_sum_utils.mod.o
	@true
g_loc_util_utils.mod:g_loc_util_utils.mod.o
	@true
glopar_utils.mod:glopar_utils.mod.o
	@true
gmopts_utils.mod:gmopts_utils.mod.o
	@true
gndstate_p_utils.mod:gndstate_p_utils.mod.o
	@true
graden_utils.mod:graden_utils.mod.o
	@true
gs_disortho_utils.mod:gs_disortho_utils.mod.o
	@true
gsize_utils.mod:gsize_utils.mod.o
	@true
gsortho_utils.mod:gsortho_utils.mod.o
	@true
gvec.mod:       gvec.mod.o
	@true
h0psi1_p_utils.mod:h0psi1_p_utils.mod.o
	@true
hardness_p_utils.mod:hardness_p_utils.mod.o
	@true
harm.mod:       harm.mod.o
	@true
header_utils.mod:header_utils.mod.o
	@true
head.mod:       head.mod.o
	@true
hesele_p_utils.mod:hesele_p_utils.mod.o
	@true
hesele_utils.mod:hesele_utils.mod.o
	@true
hess_eta_p_utils.mod:hess_eta_p_utils.mod.o
	@true
hessin_utils.mod:hessin_utils.mod.o
	@true
hessout_utils.mod:hessout_utils.mod.o
	@true
hessup_utils.mod:hessup_utils.mod.o
	@true
hfx_drivers.mod:hfx_drivers.mod.o
	@true
hfxmod.mod:     hfx.mod.o
	@true
hfx_utils.mod:  hfx_utils.mod.o
	@true
hipin_utils.mod:hipin_utils.mod.o
	@true
hip_utils.mod:  hip_utils.mod.o
	@true
hnlmat_utils.mod:hnlmat_utils.mod.o
	@true
hpsi_utils.mod: hpsi_utils.mod.o
	@true
htrstr_utils.mod:htrstr_utils.mod.o
	@true
hubbardu.mod:   hubbardu.mod.o
	@true
hubbardu_utils.mod:hubbardu_utils.mod.o
	@true
if_parallel.mod:if_parallel.mod.o
	@true
implhv.mod:     implhv.mod.o
	@true
initclust_utils.mod:initclust_utils.mod.o
	@true
initrun_driver.mod:initrun_driver.mod.o
	@true
initrun_utils.mod:initrun_utils.mod.o
	@true
inr_dr_utils.mod:inr_dr_utils.mod.o
	@true
inscan_utils.mod:inscan_utils.mod.o
	@true
interaction_manno_p_utils.mod:interaction_manno_p_utils.mod.o
	@true
interaction_p_utils.mod:interaction_p_utils.mod.o
	@true
interp3d_utils.mod:interp3d_utils.mod.o
	@true
interpt_utils.mod:interpt_utils.mod.o
	@true
ions.mod:       ions.mod.o
	@true
io_utils.mod:   io_utils.mod.o
	@true
isos.mod:       isos.mod.o
	@true
jacobi_c_utils.mod:jacobi_c_utils.mod.o
	@true
jacobi_utils.mod:jacobi_utils.mod.o
	@true
jrotation_utils.mod:jrotation_utils.mod.o
	@true
k290_2_utils.mod:k290_2_utils.mod.o
	@true
k290_utils.mod: k290_utils.mod.o
	@true
kddipo_utils.mod:kddipo_utils.mod.o
	@true
k_diis_rhofix_utils.mod:k_diis_rhofix_utils.mod.o
	@true
kdpc.mod:       kdpc.mod.o
	@true
kdp_diag_utils.mod:kdp_diag_utils.mod.o
	@true
kdp.mod:        kdp.mod.o
	@true
kdpoints_utils.mod:kdpoints_utils.mod.o
	@true
kdp_prep_utils.mod:kdp_prep_utils.mod.o
	@true
kdp_rho_utils.mod:kdp_rho_utils.mod.o
	@true
kdp_stress_kin_utils.mod:kdp_stress_kin_utils.mod.o
	@true
k_forces_driver.mod:k_forces_driver.mod.o
	@true
k_forces_utils.mod:k_forces_utils.mod.o
	@true
k_hesele_utils.mod:k_hesele_utils.mod.o
	@true
kinds.mod:      kinds.mod.o
	@true
kin_energy_utils.mod:kin_energy_utils.mod.o
	@true
k_odiis_utils.mod:k_odiis_utils.mod.o
	@true
k_pcgrad_utils.mod:k_pcgrad_utils.mod.o
	@true
kpclean_utils.mod:kpclean_utils.mod.o
	@true
kpert_potential_p_utils.mod:kpert_potential_p_utils.mod.o
	@true
kpert_util_p_utils.mod:kpert_util_p_utils.mod.o
	@true
kpnt.mod:       kpnt.mod.o
	@true
kpts.mod:       kpts.mod.o
	@true
ksdiag_utils.mod:ksdiag_utils.mod.o
	@true
ks_ener_p_utils.mod:ks_ener_p_utils.mod.o
	@true
ksmat_dist_utils.mod:ksmat_dist_utils.mod.o
	@true
ksmatmod.mod:   ksmat.mod.o
	@true
ksmat_utils.mod:ksmat_utils.mod.o
	@true
k_updwf_utils.mod:k_updwf_utils.mod.o
	@true
lanc_phon_p_utils.mod:lanc_phon_p_utils.mod.o
	@true
latgen_utils.mod:latgen_utils.mod.o
	@true
ldosmod.mod:    ldos.mod.o
	@true
ldos_utils.mod: ldos_utils.mod.o
	@true
legendre_p_utils.mod:legendre_p_utils.mod.o
	@true
linalg_utils.mod:linalg_utils.mod.o
	@true
linres.mod:     linres.mod.o
	@true
loadpa_utils.mod:loadpa_utils.mod.o
	@true
loadse_utils.mod:loadse_utils.mod.o
	@true
localize_utils.mod:localize_utils.mod.o
	@true
locpot.mod:     locpot.mod.o
	@true
lodipo_utils.mod:lodipo_utils.mod.o
	@true
lodp.mod:       lodp.mod.o
	@true
lowdin_utils.mod:lowdin_utils.mod.o
	@true
lr_diag_utils.mod:lr_diag_utils.mod.o
	@true
lr_force_utils.mod:lr_force_utils.mod.o
	@true
lr_in_utils.mod:lr_in_utils.mod.o
	@true
lr_ortho_utils.mod:lr_ortho_utils.mod.o
	@true
lr_pcg_utils.mod:lr_pcg_utils.mod.o
	@true
lr_tddft_drhoe.mod:lr_tddft_drhoe.mod.o
	@true
lr_tddft_utils.mod:lr_tddft_utils.mod.o
	@true
lr_upd_utils.mod:lr_upd_utils.mod.o
	@true
lr_xcpot_utils.mod:lr_xcpot_utils.mod.o
	@true
lscal.mod:      lscal.mod.o
	@true
lsd_elf_utils.mod:lsd_elf_utils.mod.o
	@true
lsd_func_utils.mod:lsd_func_utils.mod.o
	@true
lsfbtr_utils.mod:lsfbtr_utils.mod.o
	@true
lsforce_utils.mod:lsforce_utils.mod.o
	@true
lxc_utils.mod:  lxc_utils.mod.o
	@true
machine.mod:    machine.mod.o
	@true
matrix_p_utils.mod:matrix_p_utils.mod.o
	@true
mdclas_utils.mod:mdclas_utils.mod.o
	@true
mddiag_interaction_p_utils.mod:mddiag_interaction_p_utils.mod.o
	@true
md_driver.mod:  md_driver.mod.o
	@true
mdfile_utils.mod:mdfile_utils.mod.o
	@true
mdmain_utils.mod:mdmain_utils.mod.o
	@true
mdpt_utils.mod: mdpt_utils.mod.o
	@true
mdshop_bo_utils.mod:mdshop_bo_utils.mod.o
	@true
mdshop_cp_utils.mod:mdshop_cp_utils.mod.o
	@true
mergemod.mod:   merge.mod.o
	@true
merge_utils.mod:merge_utils.mod.o
	@true
meta_cell_utils.mod:meta_cell_utils.mod.o
	@true
meta_colvar_inp_utils.mod:meta_colvar_inp_utils.mod.o
	@true
meta_colvar_utils.mod:meta_colvar_utils.mod.o
	@true
meta_colvar_util_utils.mod:meta_colvar_util_utils.mod.o
	@true
meta_cv_qmmm_utils.mod:meta_cv_qmmm_utils.mod.o
	@true
meta_cv_utils.mod:meta_cv_utils.mod.o
	@true
meta_dyn_def_utils.mod:meta_dyn_def_utils.mod.o
	@true
meta_exlagr_methods.mod:meta_exlagr_methods.mod.o
	@true
meta_exlagr_utils.mod:meta_exlagr_utils.mod.o
	@true
meta_exl_mult_utils.mod:meta_exl_mult_utils.mod.o
	@true
meta_ex_mul_util_utils.mod:meta_ex_mul_util_utils.mod.o
	@true
metafun_utils.mod:metafun_utils.mod.o
	@true
meta_hpot_utils.mod:meta_hpot_utils.mod.o
	@true
meta_localizespin_utils.mod:meta_localizespin_utils.mod.o
	@true
meta_multiple_walkers_utils.mod:meta_multiple_walkers_utils.mod.o
	@true
metr.mod:       metr.mod.o
	@true
mfep.mod:       mfep.mod.o
	@true
min_heap.mod:   min_heap.mod.o
	@true
mixing_g_utils.mod:mixing_g_utils.mod.o
	@true
mixing_r_utils.mod:mixing_r_utils.mod.o
	@true
mltfft_utils.mod:mltfft_utils.mod.o
	@true
mm_cpmd_add_MM_forces_f77_utils.mod:mm_cpmd_add_MM_forces_f77_utils.mod.o
	@true
mm_cpmd_esp_charges_f77_utils.mod:mm_cpmd_esp_charges_f77_utils.mod.o
	@true
mm_cpmd_ext_pot_f77_utils.mod:mm_cpmd_ext_pot_f77_utils.mod.o
	@true
mm_dimmod.mod:  mm_dim.mod.o
	@true
mm_dim_utils.mod:mm_dim_utils.mod.o
	@true
mm_extrap.mod:  mm_extrap.mod.o
	@true
mm_forcematch_utils.mod:mm_forcematch_utils.mod.o
	@true
mm_forces_diag_utils.mod:mm_forces_diag_utils.mod.o
	@true
mm_forces_prop_utils.mod:mm_forces_prop_utils.mod.o
	@true
mm_init_utils.mod:mm_init_utils.mod.o
	@true
mm_input.mod:   mm_input.mod.o
	@true
mm_ion_dens.mod:mm_ion_dens.mod.o
	@true
mm_mddiag_utils.mod:mm_mddiag_utils.mod.o
	@true
mm_mdmain_utils.mod:mm_mdmain_utils.mod.o
	@true
mm_mdshop_bo_utils.mod:mm_mdshop_bo_utils.mod.o
	@true
mm_mdshop_cp_utils.mod:mm_mdshop_cp_utils.mod.o
	@true
mm_parallel.mod:mm_parallel.mod.o
	@true
mm_qmmm_forcedr_bs_utils.mod:mm_qmmm_forcedr_bs_utils.mod.o
	@true
mm_qmmm_forcedr_utils.mod:mm_qmmm_forcedr_utils.mod.o
	@true
mm_rho_forcedr_utils.mod:mm_rho_forcedr_utils.mod.o
	@true
molorb_utils.mod:molorb_utils.mod.o
	@true
mols.mod:       mols.mod.o
	@true
molstates_utils.mod:molstates_utils.mod.o
	@true
molsym_utils.mod:molsym_utils.mod.o
	@true
moverho_utils.mod:moverho_utils.mod.o
	@true
movi.mod:       movi.mod.o
	@true
mp_interface.mod:mp_interface.mod.o
	@true
mp_multiple_comm_init.mod:mp_multiple_comm_init.mod.o
	@true
mtin_utils.mod: mtin_utils.mod.o
	@true
mulliken_utils.mod:mulliken_utils.mod.o
	@true
multtb_utils.mod:multtb_utils.mod.o
	@true
mw.mod:         mw.mod.o
	@true
nabdy_ampli.mod:nabdy_ampli.mod.o
	@true
nabdy_forces.mod:nabdy_forces.mod.o
	@true
nabdy_initialize.mod:nabdy_initialize.mod.o
	@true
nabdy_md.mod:   nabdy_md.mod.o
	@true
nabdy_types.mod:nabdy_types.mod.o
	@true
newcell_utils.mod:newcell_utils.mod.o
	@true
newd_utils.mod: newd_utils.mod.o
	@true
nfunc_utils.mod:nfunc_utils.mod.o
	@true
nlcc.mod:       nlcc.mod.o
	@true
nlccset_utils.mod:nlccset_utils.mod.o
	@true
nlccstr_utils.mod:nlccstr_utils.mod.o
	@true
nlforce_utils.mod:nlforce_utils.mod.o
	@true
nlps.mod:       nlps.mod.o
	@true
nl_res_utils.mod:nl_res_utils.mod.o
	@true
nlsl_utils.mod: nlsl_utils.mod.o
	@true
nlsm1_s_utils.mod:nlsm1_s_utils.mod.o
	@true
nmr_chi_p_utils.mod:nmr_chi_p_utils.mod.o
	@true
nmr_current_p_utils.mod:nmr_current_p_utils.mod.o
	@true
nmr_full_p_utils.mod:nmr_full_p_utils.mod.o
	@true
nmr_para_p_utils.mod:nmr_para_p_utils.mod.o
	@true
nmr_position_p_utils.mod:nmr_position_p_utils.mod.o
	@true
nmr_p_utils.mod:nmr_p_utils.mod.o
	@true
nmr_shift_p_utils.mod:nmr_shift_p_utils.mod.o
	@true
nmr_util_p_utils.mod:nmr_util_p_utils.mod.o
	@true
nofo.mod:       nofo.mod.o
	@true
noforce_utils.mod:noforce_utils.mod.o
	@true
norhoe_utils.mod:norhoe_utils.mod.o
	@true
norm.mod:       norm.mod.o
	@true
nort.mod:       nort.mod.o
	@true
nosalloc_utils.mod:nosalloc_utils.mod.o
	@true
noscinit_utils.mod:noscinit_utils.mod.o
	@true
noseinit_utils.mod:noseinit_utils.mod.o
	@true
nose.mod:       nose.mod.o
	@true
noseng_utils.mod:noseng_utils.mod.o
	@true
nosepa_utils.mod:nosepa_utils.mod.o
	@true
noseup_utils.mod:noseup_utils.mod.o
	@true
nospinit_utils.mod:nospinit_utils.mod.o
	@true
npt_md_utils.mod:npt_md_utils.mod.o
	@true
nuclear_p_utils.mod:nuclear_p_utils.mod.o
	@true
numpw_utils.mod:numpw_utils.mod.o
	@true
nvarmod.mod:    nvar.mod.o
	@true
nvtx_interfaces.mod:nvtx_interfaces.mod.o
	@true
nvtx_utils.mod: nvtx_utils.mod.o
	@true
odiis_p_utils.mod:odiis_p_utils.mod.o
	@true
odiis_utils.mod:odiis_utils.mod.o
	@true
ohfd_utils.mod: ohfd_utils.mod.o
	@true
ohlr_utils.mod: ohlr_utils.mod.o
	@true
opeigr_c_utils.mod:opeigr_c_utils.mod.o
	@true
opeigr_p_utils.mod:opeigr_p_utils.mod.o
	@true
opeigr_utils.mod:opeigr_utils.mod.o
	@true
opt_lr_utils.mod:opt_lr_utils.mod.o
	@true
orbhard_utils.mod:orbhard_utils.mod.o
	@true
orbrot_utils.mod:orbrot_utils.mod.o
	@true
ortho_utils.mod:ortho_utils.mod.o
	@true
ovlap_utils.mod:ovlap_utils.mod.o
	@true
parac.mod:      parac.mod.o
	@true
para_global.mod:para_global.mod.o
	@true
part_1d.mod:    part_1d.mod.o
	@true
pbc_utils.mod:  pbc_utils.mod.o
	@true
pcgrad_driver.mod:pcgrad_driver.mod.o
	@true
pcgrad_p_utils.mod:pcgrad_p_utils.mod.o
	@true
pcgrad_utils.mod:pcgrad_utils.mod.o
	@true
pert_kpoint_p_utils.mod:pert_kpoint_p_utils.mod.o
	@true
perturbation_p_utils.mod:perturbation_p_utils.mod.o
	@true
phfac_utils.mod:phfac_utils.mod.o
	@true
phonons_p_utils.mod:phonons_p_utils.mod.o
	@true
pi_cntl_utils.mod:pi_cntl_utils.mod.o
	@true
pi_diag_utils.mod:pi_diag_utils.mod.o
	@true
pi_init_utils.mod:pi_init_utils.mod.o
	@true
pimd.mod:       pimd.mod.o
	@true
pi_mdpt_utils.mod:pi_mdpt_utils.mod.o
	@true
pi_md_utils.mod:pi_md_utils.mod.o
	@true
pimd_utils.mod: pimd_utils.mod.o
	@true
pinmtrans_utils.mod:pinmtrans_utils.mod.o
	@true
pi_npt_bomd_utils.mod:pi_npt_bomd_utils.mod.o
	@true
pi_npt_cpmd_utils.mod:pi_npt_cpmd_utils.mod.o
	@true
pi_prpt_utils.mod:pi_prpt_utils.mod.o
	@true
pi_stress_utils.mod:pi_stress_utils.mod.o
	@true
pi_wf_utils.mod:pi_wf_utils.mod.o
	@true
pm_cntl_utils.mod:pm_cntl_utils.mod.o
	@true
pm_gmopts_utils.mod:pm_gmopts_utils.mod.o
	@true
pm_init_utils.mod:pm_init_utils.mod.o
	@true
pm_mdpt_utils.mod:pm_mdpt_utils.mod.o
	@true
pm_wf_utils.mod:pm_wf_utils.mod.o
	@true
pnosmove_utils.mod:pnosmove_utils.mod.o
	@true
poin.mod:       poin.mod.o
	@true
pola.mod:       pola.mod.o
	@true
polarise_utils.mod:polarise_utils.mod.o
	@true
posupa_utils.mod:posupa_utils.mod.o
	@true
posupi_utils.mod:posupi_utils.mod.o
	@true
potfor_utils.mod:potfor_utils.mod.o
	@true
potmed_utils.mod:potmed_utils.mod.o
	@true
ppener_utils.mod:ppener_utils.mod.o
	@true
prbomd_utils.mod:prbomd_utils.mod.o
	@true
prcnosmove_utils.mod:prcnosmove_utils.mod.o
	@true
prcpmd_utils.mod:prcpmd_utils.mod.o
	@true
prcp.mod:       prcp.mod.o
	@true
prden.mod:      prden.mod.o
	@true
prep_forcematch_utils.mod:prep_forcematch_utils.mod.o
	@true
printave_utils.mod:printave_utils.mod.o
	@true
printfor_utils.mod:printfor_utils.mod.o
	@true
printpmod.mod:  printp.mod.o
	@true
printp_utils.mod:printp_utils.mod.o
	@true
prmdfile_utils.mod:prmdfile_utils.mod.o
	@true
prmem_utils.mod:prmem_utils.mod.o
	@true
prng.mod:       prng.mod.o
	@true
prng_utils.mod: prng_utils.mod.o
	@true
proja_utils.mod:proja_utils.mod.o
	@true
projv_utils.mod:projv_utils.mod.o
	@true
propin_utils.mod:propin_utils.mod.o
	@true
prop.mod:       prop.mod.o
	@true
proppt_utils.mod:proppt_utils.mod.o
	@true
prowfn_utils.mod:prowfn_utils.mod.o
	@true
proylm_utils.mod:proylm_utils.mod.o
	@true
prpcmove_utils.mod:prpcmove_utils.mod.o
	@true
prpcnosmove_utils.mod:prpcnosmove_utils.mod.o
	@true
prpnosmove_utils.mod:prpnosmove_utils.mod.o
	@true
prpt_utils.mod: prpt_utils.mod.o
	@true
prtgyr_utils.mod:prtgyr_utils.mod.o
	@true
pslo.mod:       pslo.mod.o
	@true
pstat.mod:      pstat.mod.o
	@true
ptheory_utils.mod:ptheory_utils.mod.o
	@true
purge_utils.mod:purge_utils.mod.o
	@true
putbet_utils.mod:putbet_utils.mod.o
	@true
puttau_utils.mod:puttau_utils.mod.o
	@true
pw_hfx_input_cnst.mod:pw_hfx_input_cnst.mod.o
	@true
pw_hfx.mod:     pw_hfx.mod.o
	@true
pw_hfx_resp.mod:pw_hfx_resp.mod.o
	@true
pw_hfx_resp_types.mod:pw_hfx_resp_types.mod.o
	@true
pw_hfx_resp_utils.mod:pw_hfx_resp_utils.mod.o
	@true
qrada_s_utils.mod:qrada_s_utils.mod.o
	@true
qspl.mod:       qspl.mod.o
	@true
quenbo_utils.mod:quenbo_utils.mod.o
	@true
qvan1_utils.mod:qvan1_utils.mod.o
	@true
qvan2_utils.mod:qvan2_utils.mod.o
	@true
radin_utils.mod:radin_utils.mod.o
	@true
ragg.mod:       ragg.mod.o
	@true
raman_p_utils.mod:raman_p_utils.mod.o
	@true
ranc_utils.mod: ranc_utils.mod.o
	@true
randtowf_utils.mod:randtowf_utils.mod.o
	@true
ranp_utils.mod: ranp_utils.mod.o
	@true
ratom_utils.mod:ratom_utils.mod.o
	@true
rattle_utils.mod:rattle_utils.mod.o
	@true
rbfgs_utils.mod:rbfgs_utils.mod.o
	@true
readff_utils.mod:readff_utils.mod.o
	@true
readmod.mod:    readmod.mod.o
	@true
read_prop_utils.mod:read_prop_utils.mod.o
	@true
readsr_utils.mod:readsr_utils.mod.o
	@true
readvan_utils.mod:readvan_utils.mod.o
	@true
recpnew_utils.mod:recpnew_utils.mod.o
	@true
recpupf_utils.mod:recpupf_utils.mod.o
	@true
reigs_utils.mod:reigs_utils.mod.o
	@true
rekine_utils.mod:rekine_utils.mod.o
	@true
repgen_utils.mod:repgen_utils.mod.o
	@true
resetac_utils.mod:resetac_utils.mod.o
	@true
reshaper.mod:   reshaper.mod.o
	@true
respin_p_utils.mod:respin_p_utils.mod.o
	@true
response_pmod.mod:response_p.mod.o
	@true
response_p_utils.mod:response_p_utils.mod.o
	@true
restart_p_utils.mod:restart_p_utils.mod.o
	@true
rgdiis_utils.mod:rgdiis_utils.mod.o
	@true
rggen_utils.mod:rggen_utils.mod.o
	@true
rgmopt_utils.mod:rgmopt_utils.mod.o
	@true
rgs_utils.mod:  rgs_utils.mod.o
	@true
rgsvan_utils.mod:rgsvan_utils.mod.o
	@true
rho1ofr_utils.mod:rho1ofr_utils.mod.o
	@true
rho1pri_utils.mod:rho1pri_utils.mod.o
	@true
rhodiis_utils.mod:rhodiis_utils.mod.o
	@true
rhoofr_c_utils.mod:rhoofr_c_utils.mod.o
	@true
rhoofr_kdp_utils.mod:rhoofr_kdp_utils.mod.o
	@true
rhoofr_p_utils.mod:rhoofr_p_utils.mod.o
	@true
rhoofr_utils.mod:rhoofr_utils.mod.o
	@true
rhopri_utils.mod:rhopri_utils.mod.o
	@true
rhov1_utils.mod:rhov1_utils.mod.o
	@true
rhov_utils.mod: rhov_utils.mod.o
	@true
rinforce_nuc_utils.mod:rinforce_nuc_utils.mod.o
	@true
rinforce_utils.mod:rinforce_utils.mod.o
	@true
rinit_utils.mod:rinit_utils.mod.o
	@true
rinitwf_driver.mod:rinitwf_driver.mod.o
	@true
rinitwf_utils.mod:rinitwf_utils.mod.o
	@true
rinvel_utils.mod:rinvel_utils.mod.o
	@true
rk4ov_utils.mod:rk4ov_utils.mod.o
	@true
rkpnt_utils.mod:rkpnt_utils.mod.o
	@true
rlbfgs_io.mod:  rlbfgs_io.mod.o
	@true
rlbfgs_utils.mod:rlbfgs_utils.mod.o
	@true
rmas.mod:       rmas.mod.o
	@true
rnl_dk_p_utils.mod:rnl_dk_p_utils.mod.o
	@true
rnlfl_utils.mod:rnlfl_utils.mod.o
	@true
rnlfor_utils.mod:rnlfor_utils.mod.o
	@true
rnlin_utils.mod:rnlin_utils.mod.o
	@true
rnlrh_utils.mod:rnlrh_utils.mod.o
	@true
rnlset_utils.mod:rnlset_utils.mod.o
	@true
rnlsm1_utils.mod:rnlsm1_utils.mod.o
	@true
rnlsm_2d_utils.mod:rnlsm_2d_utils.mod.o
	@true
rnlsm2_utils.mod:rnlsm2_utils.mod.o
	@true
rnlsmd_utils.mod:rnlsmd_utils.mod.o
	@true
rnlsm_p_utils.mod:rnlsm_p_utils.mod.o
	@true
rnlsm_utils.mod:rnlsm_utils.mod.o
	@true
ropt.mod:       ropt.mod.o
	@true
rortog_utils.mod:rortog_utils.mod.o
	@true
rortv_utils.mod:rortv_utils.mod.o
	@true
rotate_my_wannier_manno_p_utils.mod:rotate_my_wannier_manno_p_utils.mod.o
	@true
rotate_my_wannier_para_p_utils.mod:rotate_my_wannier_para_p_utils.mod.o
	@true
rotate_utils.mod:rotate_utils.mod.o
	@true
rotvel_utils.mod:rotvel_utils.mod.o
	@true
rpiiint_utils.mod:rpiiint_utils.mod.o
	@true
rrandd_utils.mod:rrandd_utils.mod.o
	@true
rrane_utils.mod:rrane_utils.mod.o
	@true
rreadf_utils.mod:rreadf_utils.mod.o
	@true
rrfo_utils.mod: rrfo_utils.mod.o
	@true
rscpot_utils.mod:rscpot_utils.mod.o
	@true
rscve_utils.mod:rscve_utils.mod.o
	@true
rscvp_utils.mod:rscvp_utils.mod.o
	@true
rswfmod.mod:    rswf.mod.o
	@true
rv30_utils.mod: rv30_utils.mod.o
	@true
rwfopt_nuc_utils.mod:rwfopt_nuc_utils.mod.o
	@true
rwfopt_p_utils.mod:rwfopt_p_utils.mod.o
	@true
rwfopt_utils.mod:rwfopt_utils.mod.o
	@true
rw_linres_utils.mod:rw_linres_utils.mod.o
	@true
rwswap_utils.mod:rwswap_utils.mod.o
	@true
sample_utils.mod:sample_utils.mod.o
	@true
saop_utils.mod: saop_utils.mod.o
	@true
scrp.mod:       scrp.mod.o
	@true
sdcell_utils.mod:sdcell_utils.mod.o
	@true
sd_ii_utils.mod:sd_ii_utils.mod.o
	@true
sdion_utils.mod:sdion_utils.mod.o
	@true
sdlinres_utils.mod:sdlinres_utils.mod.o
	@true
sd_loc2_utils.mod:sd_loc2_utils.mod.o
	@true
sd_loc_utils.mod:sd_loc_utils.mod.o
	@true
sd_nl2_utils.mod:sd_nl2_utils.mod.o
	@true
sd_nl_utils.mod:sd_nl_utils.mod.o
	@true
sd_wannier_utils.mod:sd_wannier_utils.mod.o
	@true
secder_driver.mod:secder_driver.mod.o
	@true
secder_utils.mod:secder_utils.mod.o
	@true
secdpt_utils.mod:secdpt_utils.mod.o
	@true
setbasis_utils.mod:setbasis_utils.mod.o
	@true
setbsstate_utils.mod:setbsstate_utils.mod.o
	@true
setcnst_utils.mod:setcnst_utils.mod.o
	@true
set_cp_grp_utils.mod:set_cp_grp_utils.mod.o
	@true
setirec_utils.mod:setirec_utils.mod.o
	@true
setsc_utils.mod:setsc_utils.mod.o
	@true
setsys_utils.mod:setsys_utils.mod.o
	@true
sfac.mod:       sfac.mod.o
	@true
sgpp.mod:       sgpp.mod.o
	@true
shake_utils.mod:shake_utils.mod.o
	@true
shock.mod:      shock.mod.o
	@true
shop_adds_utils.mod:shop_adds_utils.mod.o
	@true
shop_ekinqm.mod:shop_ekinqm.mod.o
	@true
shop.mod:       shop.mod.o
	@true
shop_rest_2.mod:shop_rest_2.mod.o
	@true
shop_rest.mod:  shop_rest.mod.o
	@true
sh_tddft_utils.mod:sh_tddft_utils.mod.o
	@true
sh_utils.mod:   sh_utils.mod.o
	@true
simple_model_p_utils.mod:simple_model_p_utils.mod.o
	@true
simulmod.mod:   simul.mod.o
	@true
sizeof_kinds.mod:sizeof_kinds.mod.o
	@true
soc.mod:        soc.mod.o
	@true
soc_types.mod:  soc_types.mod.o
	@true
softex_utils.mod:softex_utils.mod.o
	@true
soft.mod:       soft.mod.o
	@true
sort_utils.mod: sort_utils.mod.o
	@true
special_functions.mod:special_functions.mod.o
	@true
specpt_utils.mod:specpt_utils.mod.o
	@true
sphe.mod:       sphe.mod.o
	@true
spin.mod:       spin.mod.o
	@true
spsi_utils.mod: spsi_utils.mod.o
	@true
ssic_utils.mod: ssic_utils.mod.o
	@true
stagetrans_utils.mod:stagetrans_utils.mod.o
	@true
startpa_utils.mod:startpa_utils.mod.o
	@true
state_utils.mod:state_utils.mod.o
	@true
stcop_utils.mod:stcop_utils.mod.o
	@true
store_types.mod:store_types.mod.o
	@true
str2.mod:       str2.mod.o
	@true
stress_utils.mod:stress_utils.mod.o
	@true
string_utils.mod:string_utils.mod.o
	@true
strs.mod:       strs.mod.o
	@true
struc.mod:      struc.mod.o
	@true
struc_utils.mod:struc_utils.mod.o
	@true
sumfnl_utils.mod:sumfnl_utils.mod.o
	@true
summat_utils.mod:summat_utils.mod.o
	@true
swap.mod:       swap.mod.o
	@true
symm4.mod:      symm4.mod.o
	@true
symmetry_utils.mod:symmetry_utils.mod.o
	@true
symm.mod:       symm.mod.o
	@true
symtrz_utils.mod:symtrz_utils.mod.o
	@true
syscomb_utils.mod:syscomb_utils.mod.o
	@true
sysin_utils.mod:sysin_utils.mod.o
	@true
system.mod:     system.mod.o
	@true
system_utils.mod:system_utils.mod.o
	@true
tauf.mod:       tauf.mod.o
	@true
tauofr_utils.mod:tauofr_utils.mod.o
	@true
tbxc.mod:       tbxc.mod.o
	@true
td_cayley_utils.mod:td_cayley_utils.mod.o
	@true
td_dav_utils.mod:td_dav_utils.mod.o
	@true
td_force_utils.mod:td_force_utils.mod.o
	@true
td_input.mod:   td_input.mod.o
	@true
td_mm_qmmm_forcedr_utils.mod:td_mm_qmmm_forcedr_utils.mod.o
	@true
td_nacvs_utils.mod:td_nacvs_utils.mod.o
	@true
td_nhdav_utils.mod:td_nhdav_utils.mod.o
	@true
tdnlfor_utils.mod:tdnlfor_utils.mod.o
	@true
td_os_berry_utils.mod:td_os_berry_utils.mod.o
	@true
td_os_utils.mod:td_os_utils.mod.o
	@true
td_pcg_utils.mod:td_pcg_utils.mod.o
	@true
td_prop_utils.mod:td_prop_utils.mod.o
	@true
td_utils.mod:   td_utils.mod.o
	@true
temps.mod:      temps.mod.o
	@true
testex_utils.mod:testex_utils.mod.o
	@true
teststore_utils.mod:teststore_utils.mod.o
	@true
thread_view_types.mod:thread_view_types.mod.o
	@true
thread_view_utils.mod:thread_view_utils.mod.o
	@true
time.mod:       time.mod.o
	@true
timer.mod:      timer.mod.o
	@true
totstr_utils.mod:totstr_utils.mod.o
	@true
tpar.mod:       tpar.mod.o
	@true
tpot.mod:       tpot.mod.o
	@true
transmemod.mod: transme.mod.o
	@true
transme_utils.mod:transme_utils.mod.o
	@true
tst2min_inp_utils.mod:tst2min_inp_utils.mod.o
	@true
tst2min_utils.mod:tst2min_utils.mod.o
	@true
up3_p_utils.mod:up3_p_utils.mod.o
	@true
updrho_utils.mod:updrho_utils.mod.o
	@true
updwf_p_utils.mod:updwf_p_utils.mod.o
	@true
updwf_utils.mod:updwf_utils.mod.o
	@true
util_p_utils.mod:util_p_utils.mod.o
	@true
utils.mod:      utils.mod.o
	@true
u_upd_exp_sum_utils.mod:u_upd_exp_sum_utils.mod.o
	@true
u_upd_exp_utils.mod:u_upd_exp_utils.mod.o
	@true
u_upd_spread_sum_utils.mod:u_upd_spread_sum_utils.mod.o
	@true
v1ofrho1_utils.mod:v1ofrho1_utils.mod.o
	@true
v1ofrho_p_utils.mod:v1ofrho_p_utils.mod.o
	@true
v1xc_p_utils.mod:v1xc_p_utils.mod.o
	@true
vbeta_utils.mod:vbeta_utils.mod.o
	@true
vdbinit_utils.mod:vdbinit_utils.mod.o
	@true
vdbp.mod:       vdbp.mod.o
	@true
vdbt.mod:       vdbt.mod.o
	@true
vdwcmod.mod:    vdwcmod.mod.o
	@true
vdwcmod_utils.mod:vdwcmod_utils.mod.o
	@true
vdwin_utils.mod:vdwin_utils.mod.o
	@true
vdw_utils.mod:  vdw_utils.mod.o
	@true
vdw_wf_alloc_utils.mod:vdw_wf_alloc_utils.mod.o
	@true
velocitinp_utils.mod:velocitinp_utils.mod.o
	@true
velupa_utils.mod:velupa_utils.mod.o
	@true
velupi_utils.mod:velupi_utils.mod.o
	@true
vepsup_utils.mod:vepsup_utils.mod.o
	@true
vgsortho_utils.mod:vgsortho_utils.mod.o
	@true
vhk_utils.mod:  vhk_utils.mod.o
	@true
vibana_utils.mod:vibana_utils.mod.o
	@true
vlocst_utils.mod:vlocst_utils.mod.o
	@true
voa_p_utils.mod:voa_p_utils.mod.o
	@true
vofrhoa_utils.mod:vofrhoa_utils.mod.o
	@true
vofrhob_utils.mod:vofrhob_utils.mod.o
	@true
vofrhoc_utils.mod:vofrhoc_utils.mod.o
	@true
vofrhoh_utils.mod:vofrhoh_utils.mod.o
	@true
vofrhos_utils.mod:vofrhos_utils.mod.o
	@true
vofrhot_utils.mod:vofrhot_utils.mod.o
	@true
vofrho_utils.mod:vofrho_utils.mod.o
	@true
vpsi_lse_utils.mod:vpsi_lse_utils.mod.o
	@true
vpsi_p_utils.mod:vpsi_p_utils.mod.o
	@true
vpsi_utils.mod: vpsi_utils.mod.o
	@true
vtaupsi_utils.mod:vtaupsi_utils.mod.o
	@true
vtd2_utils.mod: vtd2_utils.mod.o
	@true
wannier_center_utils.mod:wannier_center_utils.mod.o
	@true
wannier_print_utils.mod:wannier_print_utils.mod.o
	@true
wann.mod:       wann.mod.o
	@true
wc_dos_utils.mod:wc_dos_utils.mod.o
	@true
wfnio_utils.mod:wfnio_utils.mod.o
	@true
wfn_print_utils.mod:wfn_print_utils.mod.o
	@true
wfopts_utils.mod:wfopts_utils.mod.o
	@true
wr30wfn_utils.mod:wr30wfn_utils.mod.o
	@true
wrccfl_utils.mod:wrccfl_utils.mod.o
	@true
wrener_utils.mod:wrener_utils.mod.o
	@true
wrgeo_utils.mod:wrgeo_utils.mod.o
	@true
wrintf_utils.mod:wrintf_utils.mod.o
	@true
write_pp_utils.mod:write_pp_utils.mod.o
	@true
wr_temps_utils.mod:wr_temps_utils.mod.o
	@true
wv30_utils.mod: wv30_utils.mod.o
	@true
xcener_utils.mod:xcener_utils.mod.o
	@true
xcstr_utils.mod:xcstr_utils.mod.o
	@true
x_hjs.mod:      x_hjs.mod.o
	@true
xinr.mod:       xinr.mod.o
	@true
ylmr2_utils.mod:ylmr2_utils.mod.o
	@true
ylmr_utils.mod: ylmr_utils.mod.o
	@true
zdiis_utils.mod:zdiis_utils.mod.o
	@true
zeroing_utils.mod:zeroing_utils.mod.o
	@true
znum_mat_utils.mod:znum_mat_utils.mod.o
	@true
################################################################################
# Object dependencies: CPMD
# 
aainit_utils.mod.f90:$(SRCDIR)/aainit_utils.mod.F90
aainit_utils.mod.o:aainit_utils.mod.f90 aavan.mod cnst.mod \
                error_handling.mod kinds.mod parac.mod system.mod \
                zeroing_utils.mod

aavan.mod.f90:  $(SRCDIR)/aavan.mod.F90
aavan.mod.o:    aavan.mod.f90 kinds.mod system.mod

adapttol_utils.mod.f90:$(SRCDIR)/adapttol_utils.mod.F90
adapttol_utils.mod.o:adapttol_utils.mod.f90 kinds.mod lscal.mod \
                mp_interface.mod parac.mod system.mod

adat.mod.f90:   $(SRCDIR)/adat.mod.F90
adat.mod.o:     adat.mod.f90 kinds.mod

adjmu_utils.mod.f90:$(SRCDIR)/adjmu_utils.mod.F90
adjmu_utils.mod.o:adjmu_utils.mod.f90 error_handling.mod fint.mod \
                kinds.mod parac.mod timer.mod

afbdr_utils.mod.f90:$(SRCDIR)/afbdr_utils.mod.F90
afbdr_utils.mod.o:afbdr_utils.mod.f90 fnlalloc_utils.mod kinds.mod \
                rnlsm_utils.mod system.mod tdnlfor_utils.mod

ainitwf_utils.mod.f90:$(SRCDIR)/ainitwf_utils.mod.F90
ainitwf_utils.mod.o:ainitwf_utils.mod.f90 atomwf_utils.mod \
                atrho_utils.mod atwf.mod error_handling.mod \
                gsortho_utils.mod kinds.mod kpts.mod ksmat_dist_utils.mod \
                ksmat_utils.mod mp_interface.mod ortho_utils.mod \
                ovlap_utils.mod parac.mod phfac_utils.mod prng_utils.mod \
                pslo.mod rscpot_utils.mod spin.mod summat_utils.mod \
                system.mod timer.mod wfnio_utils.mod zeroing_utils.mod

anderson_utils.mod.f90:$(SRCDIR)/anderson_utils.mod.F90
anderson_utils.mod.o:anderson_utils.mod.f90 andr.mod kinds.mod \
                mp_interface.mod parac.mod

andp.mod.f90:   $(SRCDIR)/andp.mod.F90
andp.mod.o:     andp.mod.f90 kinds.mod

andr.mod.f90:   $(SRCDIR)/andr.mod.F90
andr.mod.o:     andr.mod.f90 kinds.mod

anneal_utils.mod.f90:$(SRCDIR)/anneal_utils.mod.F90
anneal_utils.mod.o:anneal_utils.mod.f90 bsym.mod cnst.mod ekinpp_utils.mod \
                ions.mod kinds.mod nose.mod parac.mod pimd.mod \
                system.mod tpar.mod

array_utils.mod.f90:$(SRCDIR)/array_utils.mod.F90
array_utils.mod.o:array_utils.mod.f90 error_handling.mod kinds.mod

atimes.mod.f90: $(SRCDIR)/atimes.mod.F90
atimes.mod.o:   atimes.mod.f90 kinds.mod

atimes_utils.mod.f90:$(SRCDIR)/atimes_utils.mod.F90
atimes_utils.mod.o:atimes_utils.mod.f90 atimesmod.mod error_handling.mod \
                hpsi_utils.mod kinds.mod projv_utils.mod pslo.mod \
                spin.mod system.mod

atomc_utils.mod.f90:$(SRCDIR)/atomc_utils.mod.F90
atomc_utils.mod.o:atomc_utils.mod.f90 atwf.mod cdftmod.mod \
                cnst.mod error_handling.mod fitpack_utils.mod \
                ions.mod isos.mod kinds.mod metr.mod mp_interface.mod \
                parac.mod pbc_utils.mod qspl.mod system.mod \
                timer.mod zeroing_utils.mod

atom.mod.f90:   $(SRCDIR)/atom.mod.F90
atom.mod.o:     atom.mod.f90 kinds.mod system.mod

atoms_utils.mod.f90:$(SRCDIR)/atoms_utils.mod.F90
atoms_utils.mod.o:atoms_utils.mod.f90 adat.mod kinds.mod

atomwf_utils.mod.f90:$(SRCDIR)/atomwf_utils.mod.F90
atomwf_utils.mod.o:atomwf_utils.mod.f90 atrho_utils.mod atwf.mod \
                copot_utils.mod csmat_utils.mod elct.mod ener.mod \
                error_handling.mod fnlalloc_utils.mod func.mod \
                gs_disortho_utils.mod gsortho_utils.mod ions.mod \
                kinds.mod kpclean_utils.mod kpts.mod ksmat_dist_utils.mod \
                ksmat_utils.mod mp_interface.mod nlcc.mod nlps.mod \
                ovlap_utils.mod parac.mod phfac_utils.mod prmem_utils.mod \
                pslo.mod randtowf_utils.mod reshaper.mod rgs_utils.mod \
                rnlsm_utils.mod sfac.mod sphe.mod spin.mod \
                summat_utils.mod system.mod timer.mod utils.mod \
                vofrho_utils.mod zeroing_utils.mod

atrho_utils.mod.f90:$(SRCDIR)/atrho_utils.mod.F90
atrho_utils.mod.o:atrho_utils.mod.f90 atwf.mod bessm_utils.mod \
                cppt.mod dotp_utils.mod error_handling.mod \
                fftmain_utils.mod fftnew_utils.mod fitpack_utils.mod \
                geq0mod.mod ions.mod kinds.mod mp_interface.mod \
                parac.mod qspl.mod reshaper.mod setbasis_utils.mod \
                sfac.mod spin.mod system.mod timer.mod zeroing_utils.mod

atwf.mod.f90:   $(SRCDIR)/atwf.mod.F90
atwf.mod.o:     atwf.mod.f90 kinds.mod system.mod

augchg_utils.mod.f90:$(SRCDIR)/augchg_utils.mod.F90
augchg_utils.mod.o:augchg_utils.mod.f90 cvan.mod error_handling.mod \
                ions.mod kinds.mod nlps.mod pslo.mod system.mod

azzero_utils.mod.f90:$(SRCDIR)/azzero_utils.mod.F90
azzero_utils.mod.o:azzero_utils.mod.f90 kinds.mod

bc.mod.f90:     $(SRCDIR)/bc.mod.F90
bc.mod.o:       bc.mod.f90 kinds.mod

benc.mod.f90:   $(SRCDIR)/benc.mod.F90
benc.mod.o:     benc.mod.f90

bessm_utils.mod.f90:$(SRCDIR)/bessm_utils.mod.F90
bessm_utils.mod.o:bessm_utils.mod.f90 error_handling.mod kinds.mod \
                parac.mod

bicanonicalCalculationConfig.mod.f90:$(SRCDIR)/bicanonicalCalculationConfig.mod.F90
bicanonicalCalculationConfig.mod.o:bicanonicalCalculationConfig.mod.f90 \
                kinds.mod

bicanonicalCalculation.mod.f90:$(SRCDIR)/bicanonicalCalculation.mod.F90
bicanonicalCalculation.mod.o:bicanonicalCalculation.mod.f90 \
                kinds.mod bicanonicalCalculationConfig.mod

bicanonicalConfig.mod.f90:$(SRCDIR)/bicanonicalConfig.mod.F90
bicanonicalConfig.mod.o:bicanonicalConfig.mod.f90 error_handling.mod \
                bicanonicalCalculationConfig.mod

bicanonicalCpmd.mod.f90:$(SRCDIR)/bicanonicalCpmd.mod.F90
bicanonicalCpmd.mod.o:bicanonicalCpmd.mod.f90 cnst.mod coor.mod \
                clas.mod adat.mod dum2_utils.mod elct.mod error_handling.mod \
                fillc_utils.mod fileopen_utils.mod fileopenmod.mod \
                filnmod.mod ions.mod kinds.mod mm_dimmod.mod \
                mp_interface.mod parac.mod pbc_utils.mod pimd.mod \
                readsr_utils.mod set_cp_grp_utils.mod soft.mod \
                system.mod testex_utils.mod timer.mod bicanonicalConfig.mod \
                bicanonicalCalculationConfig.mod bicanonicalCalculation.mod \
                wann.mod zeroing_utils.mod

bicanonicalInputReader.mod.f90:$(SRCDIR)/bicanonicalInputReader.mod.F90
bicanonicalInputReader.mod.o:bicanonicalInputReader.mod.f90 \
                error_handling.mod kinds.mod bicanonicalConfig.mod

blas_tuned_ES.f90:$(SRCDIR)/blas_tuned_ES.F90
blas_tuned_ES.o:blas_tuned_ES.f90 kinds.mod error_handling.mod \
                timer.mod kinds.mod error_handling.mod timer.mod \
                kinds.mod error_handling.mod timer.mod kinds.mod \
                error_handling.mod timer.mod kinds.mod error_handling.mod \
                timer.mod kinds.mod error_handling.mod timer.mod \
                kinds.mod error_handling.mod timer.mod kinds.mod \
                error_handling.mod timer.mod kinds.mod error_handling.mod \
                timer.mod kinds.mod error_handling.mod timer.mod

blas_tuned_NECSX.f90:$(SRCDIR)/blas_tuned_NECSX.F90
blas_tuned_NECSX.o:blas_tuned_NECSX.f90 kinds.mod error_handling.mod \
                timer.mod kinds.mod error_handling.mod timer.mod \
                kinds.mod error_handling.mod timer.mod kinds.mod \
                error_handling.mod timer.mod kinds.mod error_handling.mod \
                timer.mod

blas_tuned_SR11K.f90:$(SRCDIR)/blas_tuned_SR11K.F90
blas_tuned_SR11K.o:blas_tuned_SR11K.f90 kinds.mod error_handling.mod \
                timer.mod kinds.mod error_handling.mod timer.mod \
                kinds.mod error_handling.mod timer.mod kinds.mod \
                error_handling.mod timer.mod kinds.mod error_handling.mod \
                timer.mod kinds.mod error_handling.mod timer.mod \
                kinds.mod error_handling.mod timer.mod kinds.mod \
                error_handling.mod timer.mod kinds.mod error_handling.mod \
                timer.mod kinds.mod error_handling.mod timer.mod \
                kinds.mod error_handling.mod timer.mod

blas_tuned_X1.f90:$(SRCDIR)/blas_tuned_X1.F90
blas_tuned_X1.o:blas_tuned_X1.f90 kinds.mod error_handling.mod \
                timer.mod kinds.mod error_handling.mod timer.mod \
                kinds.mod error_handling.mod timer.mod

bogol_utils.mod.f90:$(SRCDIR)/bogol_utils.mod.F90
bogol_utils.mod.o:bogol_utils.mod.f90 dotp_utils.mod error_handling.mod \
                hpsi_utils.mod kinds.mod kpnt.mod kpts.mod \
                mp_interface.mod parac.mod spin.mod system.mod \
                tauf.mod

box_boundary_utils.mod.f90:$(SRCDIR)/box_boundary_utils.mod.F90
box_boundary_utils.mod.o:box_boundary_utils.mod.f90 cell.mod \
                comvel_utils.mod comvelmod.mod ions.mod isos.mod \
                kinds.mod

broyden_utils.mod.f90:$(SRCDIR)/broyden_utils.mod.F90
broyden_utils.mod.o:broyden_utils.mod.f90 broy.mod dotp_utils.mod \
                error_handling.mod geq0mod.mod kinds.mod mp_interface.mod \
                ovlap_utils.mod parac.mod spin.mod summat_utils.mod \
                timer.mod utils.mod zeroing_utils.mod

broy.mod.f90:   $(SRCDIR)/broy.mod.F90
broy.mod.o:     broy.mod.f90 kinds.mod

bs_forces_diag_utils.mod.f90:$(SRCDIR)/bs_forces_diag_utils.mod.F90
bs_forces_diag_utils.mod.o:bs_forces_diag_utils.mod.f90 bsym.mod \
                bsympnt.mod cnstfc_utils.mod cotr.mod elct.mod \
                ener.mod error_handling.mod kinds.mod lsforce_utils.mod \
                machine.mod mp_interface.mod norm.mod parac.mod \
                ropt.mod setbsstate_utils.mod soft.mod spin.mod \
                store_types.mod system.mod testex_utils.mod \
                tpar.mod updwf_utils.mod wrccfl_utils.mod wrener_utils.mod \
                wv30_utils.mod

bswfo_utils.mod.f90:$(SRCDIR)/bswfo_utils.mod.F90
bswfo_utils.mod.o:bswfo_utils.mod.f90 bsym.mod bsympnt.mod \
                coor.mod elct.mod ener.mod error_handling.mod \
                finalp_utils.mod geofile_utils.mod gsize_utils.mod \
                kinds.mod lsforce_utils.mod mm_dim_utils.mod \
                mm_dimmod.mod norm.mod parac.mod ropt.mod rwfopt_utils.mod \
                setbsstate_utils.mod setirec_utils.mod soft.mod \
                store_types.mod system.mod wrgeo_utils.mod \
                wv30_utils.mod zeroing_utils.mod

bsym.mod.f90:   $(SRCDIR)/bsym.mod.F90
bsym.mod.o:     bsym.mod.f90 kinds.mod

bsympnt.mod.f90:$(SRCDIR)/bsympnt.mod.F90
bsympnt.mod.o:  bsympnt.mod.f90 kinds.mod

calc_alm_utils.mod.f90:$(SRCDIR)/calc_alm_utils.mod.F90
calc_alm_utils.mod.o:calc_alm_utils.mod.f90 cppt.mod elct.mod \
                error_handling.mod fint.mod geq0mod.mod ions.mod \
                kinds.mod kpnt.mod kpts.mod nlps.mod ovlap_utils.mod \
                parac.mod prmem_utils.mod summat_utils.mod \
                system.mod timer.mod utils.mod

calc_pij_utils.mod.f90:$(SRCDIR)/calc_pij_utils.mod.F90
calc_pij_utils.mod.o:calc_pij_utils.mod.f90 cppt.mod geq0mod.mod \
                kinds.mod kpnt.mod kpts.mod system.mod

canon_utils.mod.f90:$(SRCDIR)/canon_utils.mod.F90
canon_utils.mod.o:canon_utils.mod.f90 error_handling.mod kinds.mod \
                mp_interface.mod ovlap_utils.mod parac.mod \
                rotate_utils.mod spin.mod summat_utils.mod \
                system.mod

cdft.mod.f90:   $(SRCDIR)/cdft.mod.F90
cdft.mod.o:     cdft.mod.f90 kinds.mod system.mod

cdft_utils.mod.f90:$(SRCDIR)/cdft_utils.mod.F90
cdft_utils.mod.o:cdft_utils.mod.f90 adat.mod atomc_utils.mod \
                atwf.mod cdftmod.mod cnst.mod coor.mod cppt.mod \
                elct.mod ener.mod error_handling.mod fft_maxfft.mod \
                fftmain_utils.mod fitpack_utils.mod ions.mod \
                isos.mod kinds.mod metr.mod mp_interface.mod \
                parac.mod pbc_utils.mod pslo.mod qspl.mod rhoofr_utils.mod \
                rnlsm_utils.mod spin.mod system.mod timer.mod \
                zeroing_utils.mod

cell.mod.f90:   $(SRCDIR)/cell.mod.F90
cell.mod.o:     cell.mod.f90 kinds.mod

chain_dr_utils.mod.f90:$(SRCDIR)/chain_dr_utils.mod.F90
chain_dr_utils.mod.o:chain_dr_utils.mod.f90 cnst_dyn.mod error_handling.mod \
                ions.mod kinds.mod metr.mod mm_dimmod.mod mm_input.mod \
                parac.mod timer.mod zeroing_utils.mod

chksym_utils.mod.f90:$(SRCDIR)/chksym_utils.mod.F90
chksym_utils.mod.o:chksym_utils.mod.f90 error_handling.mod \
                ions.mod k290_2_utils.mod kinds.mod metr.mod \
                molsym_utils.mod multtb_utils.mod parac.mod \
                readsr_utils.mod rmas.mod symm.mod symtrz_utils.mod \
                system.mod timer.mod utils.mod zeroing_utils.mod

clas_force_utils.mod.f90:$(SRCDIR)/clas_force_utils.mod.F90
clas_force_utils.mod.o:clas_force_utils.mod.f90 clas.mod error_handling.mod \
                kinds.mod mp_interface.mod parac.mod pbc_utils.mod \
                timer.mod utils.mod zeroing_utils.mod

clas.mod.f90:   $(SRCDIR)/clas.mod.F90
clas.mod.o:     clas.mod.f90 kinds.mod

clinbcg_utils.mod.f90:$(SRCDIR)/clinbcg_utils.mod.F90
clinbcg_utils.mod.o:clinbcg_utils.mod.f90 atimes_utils.mod \
                atimesmod.mod cppt.mod dotp_utils.mod error_handling.mod \
                kinds.mod kpts.mod machine.mod parac.mod system.mod

cl_init_utils.mod.f90:$(SRCDIR)/cl_init_utils.mod.F90
cl_init_utils.mod.o:cl_init_utils.mod.f90 clas.mod error_handling.mod \
                inscan_utils.mod kinds.mod mp_interface.mod \
                parac.mod readff_utils.mod readsr_utils.mod \
                rmas.mod

cmaos_utils.mod.f90:$(SRCDIR)/cmaos_utils.mod.F90
cmaos_utils.mod.o:cmaos_utils.mod.f90 atwf.mod error_handling.mod \
                ions.mod kinds.mod parac.mod prop.mod utils.mod \
                zeroing_utils.mod

c_mem_utils.o:  $(SRCDIR)/c_mem_utils.c

cnst_dyn.mod.f90:$(SRCDIR)/cnst_dyn.mod.F90
cnst_dyn.mod.o: cnst_dyn.mod.f90 kinds.mod

cnstfc_utils.mod.f90:$(SRCDIR)/cnstfc_utils.mod.F90
cnstfc_utils.mod.o:cnstfc_utils.mod.f90 constr_utils.mod cotr.mod \
                dum2_utils.mod ener.mod error_handling.mod \
                fillc_utils.mod ions.mod isos.mod kinds.mod \
                meta_cv_utils.mod parac.mod pbc_utils.mod puttau_utils.mod \
                rmas.mod system.mod timer.mod zeroing_utils.mod

cnst.mod.f90:   $(SRCDIR)/cnst.mod.F90
cnst.mod.o:     cnst.mod.f90 kinds.mod

cnstpr_utils.mod.f90:$(SRCDIR)/cnstpr_utils.mod.F90
cnstpr_utils.mod.o:cnstpr_utils.mod.f90 coninp_utils.mod cotr.mod \
                kinds.mod mm_dim_utils.mod mm_dimmod.mod parac.mod

cofor_utils.mod.f90:$(SRCDIR)/cofor_utils.mod.F90
cofor_utils.mod.o:cofor_utils.mod.f90 cppt.mod ions.mod kinds.mod \
                nlcc.mod sfac.mod system.mod

compress.f90:   $(SRCDIR)/compress.F90
compress.o:     compress.f90 kinds.mod error_handling.mod timer.mod \
                reshaper.mod kinds.mod error_handling.mod timer.mod \
                reshaper.mod string_utils.mod parac.mod string_utils.mod \
                parac.mod kinds.mod error_handling.mod timer.mod \
                kinds.mod error_handling.mod timer.mod kinds.mod \
                error_handling.mod timer.mod kinds.mod error_handling.mod \
                timer.mod

comvel.mod.f90: $(SRCDIR)/comvel.mod.F90
comvel.mod.o:   comvel.mod.f90

comvel_utils.mod.f90:$(SRCDIR)/comvel_utils.mod.F90
comvel_utils.mod.o:comvel_utils.mod.f90 cnst.mod cotr.mod ekinpp_utils.mod \
                ions.mod kinds.mod nose.mod puttau_utils.mod \
                rmas.mod system.mod

conduct_utils.mod.f90:$(SRCDIR)/conduct_utils.mod.F90
conduct_utils.mod.o:conduct_utils.mod.f90 cnst.mod condu.mod \
                cppt.mod error_handling.mod fileopen_utils.mod \
                fileopenmod.mod kinds.mod kpnt.mod kpts.mod \
                mp_interface.mod parac.mod system.mod zeroing_utils.mod

condu.mod.f90:  $(SRCDIR)/condu.mod.F90
condu.mod.o:    condu.mod.f90 kinds.mod

coninp_utils.mod.f90:$(SRCDIR)/coninp_utils.mod.F90
coninp_utils.mod.o:coninp_utils.mod.f90 cnst.mod cotr.mod error_handling.mod \
                fileopen_utils.mod fileopenmod.mod ions.mod \
                kinds.mod mm_dim_utils.mod mm_dimmod.mod mm_input.mod \
                parac.mod readsr_utils.mod zeroing_utils.mod

constr_utils.mod.f90:$(SRCDIR)/constr_utils.mod.F90
constr_utils.mod.o:constr_utils.mod.f90 kinds.mod pbc_utils.mod

control_bcast_utils.mod.f90:$(SRCDIR)/control_bcast_utils.mod.F90
control_bcast_utils.mod.o:control_bcast_utils.mod.f90 andr.mod \
                benc.mod broy.mod cdftmod.mod comvelmod.mod \
                cotr.mod cp_cuda_types.mod error_handling.mod \
                fileopenmod.mod g_loc.mod glemod.mod ions.mod \
                linres.mod mergemod.mod mm_input.mod mp_interface.mod \
                mw.mod nabdy_types.mod nort.mod nose.mod para_global.mod \
                parac.mod prden.mod qspl.mod shop.mod shop_rest.mod \
                spin.mod store_types.mod system.mod time.mod \
                vdwcmod.mod wann.mod xinr.mod

control_def_utils.mod.f90:$(SRCDIR)/control_def_utils.mod.F90
control_def_utils.mod.o:control_def_utils.mod.f90 andr.mod \
                atwf.mod benc.mod broy.mod cdftmod.mod cnst_dyn.mod \
                comvelmod.mod conv.mod cotr.mod cp_cuda_types.mod \
                fileopenmod.mod fint.mod g_loc.mod glemod.mod \
                ions.mod isos.mod kpts.mod linres.mod machine.mod \
                nort.mod nose.mod parac.mod prden.mod prop.mod \
                qspl.mod readmod.mod spin.mod store_types.mod \
                struc.mod system.mod time.mod vdwcmod.mod wann.mod \
                xinr.mod zeroing_utils.mod

control_pri_utils.mod.f90:$(SRCDIR)/control_pri_utils.mod.F90
control_pri_utils.mod.o:control_pri_utils.mod.f90 andr.mod \
                atwf.mod broy.mod comvelmod.mod cotr.mod envj.mod \
                fileopenmod.mod fint.mod glemod.mod isos.mod \
                kpts.mod lscal.mod nort.mod nose.mod parac.mod \
                qspl.mod readsr_utils.mod shop.mod store_types.mod \
                system.mod vdwcmod.mod wann.mod xinr.mod

control_test_utils.mod.f90:$(SRCDIR)/control_test_utils.mod.F90
control_test_utils.mod.o:control_test_utils.mod.f90 broy.mod \
                cdftmod.mod comvelmod.mod cp_cuda_types.mod \
                error_handling.mod fileopenmod.mod fint.mod \
                kpts.mod linres.mod lscal.mod mm_input.mod \
                mw.mod parac.mod qspl.mod spin.mod store_types.mod \
                system.mod vdwcmod.mod wann.mod

control_utils.mod.f90:$(SRCDIR)/control_utils.mod.F90
control_utils.mod.o:control_utils.mod.f90 andr.mod atwf.mod \
                benc.mod broy.mod cdftmod.mod comvelmod.mod \
                control_bcast_utils.mod control_def_utils.mod \
                control_pri_utils.mod control_test_utils.mod \
                cotr.mod cp_cuda_types.mod envj.mod error_handling.mod \
                fileopenmod.mod fint.mod g_loc.mod glemod.mod \
                header_utils.mod hubbardu.mod inscan_utils.mod \
                ions.mod isos.mod kinds.mod kpts.mod linres.mod \
                lscal.mod machine.mod mergemod.mod mm_input.mod \
                mm_parallel.mod nabdy_types.mod nort.mod nose.mod \
                para_global.mod parac.mod prden.mod qspl.mod \
                readsr_utils.mod rlbfgs_utils.mod shop.mod \
                shop_rest.mod spin.mod store_types.mod struc.mod \
                system.mod time.mod bicanonicalInputReader.mod \
                bicanonicalCpmd.mod vdwcmod.mod wann.mod xinr.mod

conv.mod.f90:   $(SRCDIR)/conv.mod.F90
conv.mod.o:     conv.mod.f90

coor.mod.f90:   $(SRCDIR)/coor.mod.F90
coor.mod.o:     coor.mod.f90 kinds.mod

copot_utils.mod.f90:$(SRCDIR)/copot_utils.mod.F90
copot_utils.mod.o:copot_utils.mod.f90 corec_utils.mod cppt.mod \
                dotp_utils.mod error_handling.mod fftmain_utils.mod \
                gcener_utils.mod graden_utils.mod ions.mod \
                kinds.mod nlcc.mod parac.mod sfac.mod spin.mod \
                strs.mod system.mod timer.mod utils.mod xcener_utils.mod \
                zeroing_utils.mod

corec_utils.mod.f90:$(SRCDIR)/corec_utils.mod.F90
corec_utils.mod.o:corec_utils.mod.f90 cppt.mod fft_maxfft.mod \
                fftmain_utils.mod fftnew_utils.mod geq0mod.mod \
                ions.mod kinds.mod nlcc.mod parac.mod sfac.mod \
                spin.mod system.mod timer.mod zeroing_utils.mod

cores.mod.f90:  $(SRCDIR)/cores.mod.F90
cores.mod.o:    cores.mod.f90 kinds.mod

core_spect_utils.mod.f90:$(SRCDIR)/core_spect_utils.mod.F90
core_spect_utils.mod.o:core_spect_utils.mod.f90 adat.mod atwf.mod \
                cnst.mod cores.mod cppt.mod error_handling.mod \
                fileopen_utils.mod fileopenmod.mod fitpack_utils.mod \
                gvec.mod ions.mod kinds.mod kpnt.mod kpts.mod \
                lsfbtr_utils.mod mp_interface.mod parac.mod \
                qspl.mod recpnew_utils.mod setbasis_utils.mod \
                sfac.mod sphe.mod system.mod ylmr_utils.mod \
                zeroing_utils.mod

cotr.mod.f90:   $(SRCDIR)/cotr.mod.F90
cotr.mod.o:     cotr.mod.f90 kinds.mod

cp_cuda_types.mod.f90:$(SRCDIR)/cp_cuda_types.mod.F90
cp_cuda_types.mod.o:cp_cuda_types.mod.f90

cp_cuda_utils.mod.f90:$(SRCDIR)/cp_cuda_utils.mod.F90
cp_cuda_utils.mod.o:cp_cuda_utils.mod.f90 cp_cuda_types.mod \
                cublas_types.mod cublas_utils.mod cuda_utils.mod \
                cufft_utils.mod error_handling.mod mp_interface.mod \
                parac.mod timer.mod

cp_cudensity_utils.mod.f90:$(SRCDIR)/cp_cudensity_utils.mod.F90
cp_cudensity_utils.mod.o:cp_cudensity_utils.mod.f90 cp_cufft_types.mod \
                cp_curho_types.mod cuda_types.mod cuda_utils.mod \
                cuuser_utils.mod kinds.mod thread_view_types.mod

cp_cufft_types.mod.f90:$(SRCDIR)/cp_cufft_types.mod.F90
cp_cufft_types.mod.o:cp_cufft_types.mod.f90 cublas_types.mod \
                cuda_types.mod cufft_types.mod error_handling.mod

cp_cufft_utils.mod.f90:$(SRCDIR)/cp_cufft_utils.mod.F90
cp_cufft_utils.mod.o:cp_cufft_utils.mod.f90 cp_cuda_types.mod \
                cp_cufft_types.mod cublas_utils.mod cuda_types.mod \
                cuda_utils.mod cufft_interfaces.mod cufft_types.mod \
                cufft_utils.mod error_handling.mod fft.mod \
                fft_maxfft.mod isos.mod kinds.mod parac.mod \
                sizeof_kinds.mod timer.mod cp_cuda_types.mod

cp_cuortho_types.mod.f90:$(SRCDIR)/cp_cuortho_types.mod.F90
cp_cuortho_types.mod.o:cp_cuortho_types.mod.f90 cublas_types.mod \
                cuda_types.mod cusolver_types.mod error_handling.mod

cp_cuortho_utils.mod.f90:$(SRCDIR)/cp_cuortho_utils.mod.F90
cp_cuortho_utils.mod.o:cp_cuortho_utils.mod.f90 cp_cuda_types.mod \
                cp_cuortho_types.mod cublas_utils.mod cuda_utils.mod \
                cusolver_utils.mod error_handling.mod kinds.mod \
                sizeof_kinds.mod timer.mod

cp_curho_types.mod.f90:$(SRCDIR)/cp_curho_types.mod.F90
cp_curho_types.mod.o:cp_curho_types.mod.f90 cuda_types.mod \
                error_handling.mod

cp_curho_utils.mod.f90:$(SRCDIR)/cp_curho_utils.mod.F90
cp_curho_utils.mod.o:cp_curho_utils.mod.f90 cp_cuda_types.mod \
                cp_curho_types.mod cuda_utils.mod error_handling.mod \
                kinds.mod sizeof_kinds.mod timer.mod

cp_cuvpsi_types.mod.f90:$(SRCDIR)/cp_cuvpsi_types.mod.F90
cp_cuvpsi_types.mod.o:cp_cuvpsi_types.mod.f90 cuda_types.mod

cp_cuvpsi_utils.mod.f90:$(SRCDIR)/cp_cuvpsi_utils.mod.F90
cp_cuvpsi_utils.mod.o:cp_cuvpsi_utils.mod.f90 cp_cuda_types.mod \
                cp_cufft_types.mod cp_cuvpsi_types.mod cuda_types.mod \
                cuda_utils.mod cuuser_utils.mod error_handling.mod \
                kinds.mod sizeof_kinds.mod thread_view_types.mod \
                timer.mod

cp_cuwfn_types.mod.f90:$(SRCDIR)/cp_cuwfn_types.mod.F90
cp_cuwfn_types.mod.o:cp_cuwfn_types.mod.f90 cuda_types.mod \
                error_handling.mod kinds.mod

cp_cuwfn_utils.mod.f90:$(SRCDIR)/cp_cuwfn_utils.mod.F90
cp_cuwfn_utils.mod.o:cp_cuwfn_utils.mod.f90 cp_cuda_types.mod \
                cp_cuwfn_types.mod cuda_utils.mod error_handling.mod \
                kinds.mod sizeof_kinds.mod timer.mod

cp_dgga_correlation_utils.mod.f90:$(SRCDIR)/cp_dgga_correlation_utils.mod.F90
cp_dgga_correlation_utils.mod.o:cp_dgga_correlation_utils.mod.f90 \
                cpfunc_types.mod func.mod kinds.mod system.mod

cp_dgga_exchange_utils.mod.f90:$(SRCDIR)/cp_dgga_exchange_utils.mod.F90
cp_dgga_exchange_utils.mod.o:cp_dgga_exchange_utils.mod.f90 \
                cpfunc_types.mod func.mod kinds.mod system.mod

cp_dxc_driver.mod.f90:$(SRCDIR)/cp_dxc_driver.mod.F90
cp_dxc_driver.mod.o:cp_dxc_driver.mod.f90 cp_xc_utils.mod cpfunc_types.mod \
                error_handling.mod kinds.mod parac.mod system.mod

cpfunc_types.mod.f90:$(SRCDIR)/cpfunc_types.mod.F90
cpfunc_types.mod.o:cpfunc_types.mod.f90 kinds.mod

cpfunc_utils.mod.f90:$(SRCDIR)/cpfunc_utils.mod.F90
cpfunc_utils.mod.o:cpfunc_utils.mod.f90

cp_gga_correlation_utils.mod.f90:$(SRCDIR)/cp_gga_correlation_utils.mod.F90
cp_gga_correlation_utils.mod.o:cp_gga_correlation_utils.mod.f90 \
                cnst.mod cp_lda_correlation_utils.mod cpfunc_types.mod \
                kinds.mod

cp_gga_exchange_utils.mod.f90:$(SRCDIR)/cp_gga_exchange_utils.mod.F90
cp_gga_exchange_utils.mod.o:cp_gga_exchange_utils.mod.f90 cnst.mod \
                cpfunc_types.mod func.mod kinds.mod

cp_grp_utils.mod.f90:$(SRCDIR)/cp_grp_utils.mod.F90
cp_grp_utils.mod.o:cp_grp_utils.mod.f90 error_handling.mod \
                geq0mod.mod kinds.mod kpts.mod mp_interface.mod \
                parac.mod system.mod timer.mod kinds.mod error_handling.mod \
                timer.mod kinds.mod error_handling.mod timer.mod \
                parac.mod kinds.mod error_handling.mod timer.mod \
                kinds.mod error_handling.mod timer.mod kinds.mod \
                error_handling.mod timer.mod

cp_ieee_interface.mod.f90:$(SRCDIR)/cp_ieee_interface.mod.F90
cp_ieee_interface.mod.o:cp_ieee_interface.mod.f90 kinds.mod

cp_lda_correlation_utils.mod.f90:$(SRCDIR)/cp_lda_correlation_utils.mod.F90
cp_lda_correlation_utils.mod.o:cp_lda_correlation_utils.mod.f90 \
                cnst.mod cpfunc_types.mod kinds.mod

cp_lda_exchange_utils.mod.f90:$(SRCDIR)/cp_lda_exchange_utils.mod.F90
cp_lda_exchange_utils.mod.o:cp_lda_exchange_utils.mod.f90 cpfunc_types.mod \
                func.mod kinds.mod

cplngs.mod.f90: $(SRCDIR)/cplngs.mod.F90
cplngs.mod.o:   cplngs.mod.f90 kinds.mod

cplngs_utils.mod.f90:$(SRCDIR)/cplngs_utils.mod.F90
cplngs_utils.mod.o:cplngs_utils.mod.f90 andp.mod canon_utils.mod \
                coor.mod cplngsmod.mod cppt.mod eicalc_utils.mod \
                eind_ii_utils.mod eind_loc_utils.mod eind_nl_utils.mod \
                elct.mod ener.mod error_handling.mod fftmain_utils.mod \
                fnlalloc_utils.mod forcedr_driver.mod forcep_utils.mod \
                forces_diag_utils.mod geq0mod.mod hpsi_utils.mod \
                inscan_utils.mod ions.mod kinds.mod kpts.mod \
                ksdiag_utils.mod linres.mod lr_xcpot_utils.mod \
                lscal.mod mp_interface.mod nl_res_utils.mod \
                nlps.mod opt_lr_utils.mod ovlap_utils.mod parac.mod \
                phfac_utils.mod poin.mod pslo.mod rho1ofr_utils.mod \
                rhoofr_utils.mod rnlsm1_utils.mod rnlsm2_utils.mod \
                rnlsm_2d_utils.mod rnlsm_utils.mod ropt.mod \
                sfac.mod sgpp.mod spin.mod summat_utils.mod \
                symtrz_utils.mod system.mod timer.mod v1ofrho1_utils.mod \
                zeroing_utils.mod

cpmd.f90:       $(SRCDIR)/cpmd.F90
cpmd.o:         cpmd.f90 kinds.mod error_handling.mod timer.mod \
                machine.mod mp_interface.mod prng_utils.mod \
                control_utils.mod dftin_utils.mod sysin_utils.mod \
                setsc_utils.mod detsp_utils.mod mm_init_utils.mod \
                read_prop_utils.mod ratom_utils.mod vdwin_utils.mod \
                propin_utils.mod respin_p_utils.mod setsys_utils.mod \
                genxc_utils.mod numpw_utils.mod pi_cntl_utils.mod \
                pi_init_utils.mod meta_multiple_walkers_utils.mod \
                bicanonicalCpmd.mod bicanonicalConfig.mod nmr_para_p_utils.mod \
                rinit_utils.mod rinforce_utils.mod fftprp_utils.mod \
                initclust_utils.mod dginit_utils.mod nosalloc_utils.mod \
                exterp_utils.mod setbasis_utils.mod dqgalloc_utils.mod \
                gle_utils.mod pi_wf_utils.mod pi_mdpt_utils.mod \
                pi_prpt_utils.mod pm_wf_utils.mod pm_gmopts_utils.mod \
                pm_mdpt_utils.mod prpt_utils.mod mdpt_utils.mod \
                gmopts_utils.mod vdw_wf_alloc_utils.mod wfopts_utils.mod \
                interpt_utils.mod secdpt_utils.mod proppt_utils.mod \
                orbhard_utils.mod response_p_utils.mod prep_forcematch_utils.mod \
                prmem_utils.mod system.mod parac.mod pimd.mod \
                mw.mod specpt_utils.mod ropt.mod soc.mod set_cp_grp_utils.mod \
                pm_init_utils.mod startpa_utils.mod envir_utils.mod \
                pm_cntl_utils.mod setcnst_utils.mod softex_utils.mod \
                fileopen_utils.mod linres.mod lr_in_utils.mod \
                cp_cuda_utils.mod ortho_utils.mod

cp_mgga_correlation_utils.mod.f90:$(SRCDIR)/cp_mgga_correlation_utils.mod.F90
cp_mgga_correlation_utils.mod.o:cp_mgga_correlation_utils.mod.f90 \
                kinds.mod cnst.mod func.mod cpfunc_types.mod \
                cp_lda_correlation_utils.mod cp_gga_correlation_utils.mod \
                cp_mgga_exchange_utils.mod

cp_mgga_exchange_utils.mod.f90:$(SRCDIR)/cp_mgga_exchange_utils.mod.F90
cp_mgga_exchange_utils.mod.o:cp_mgga_exchange_utils.mod.f90 \
                cnst.mod cpfunc_types.mod func.mod kinds.mod

cppt.mod.f90:   $(SRCDIR)/cppt.mod.F90
cppt.mod.o:     cppt.mod.f90 kinds.mod system.mod

cp_xc_driver.mod.f90:$(SRCDIR)/cp_xc_driver.mod.F90
cp_xc_driver.mod.o:cp_xc_driver.mod.f90 cp_xc_utils.mod cpfunc_types.mod \
                error_handling.mod kinds.mod lxc_utils.mod \
                parac.mod system.mod timer.mod

cp_xc_utils.mod.f90:$(SRCDIR)/cp_xc_utils.mod.F90
cp_xc_utils.mod.o:cp_xc_utils.mod.f90 cp_gga_correlation_utils.mod \
                cp_gga_exchange_utils.mod cp_dgga_exchange_utils.mod \
                cp_dgga_correlation_utils.mod cp_lda_correlation_utils.mod \
                cp_lda_exchange_utils.mod cp_mgga_exchange_utils.mod \
                cp_mgga_correlation_utils.mod cpfunc_types.mod \
                error_handling.mod func.mod kinds.mod lxc_utils.mod

crotwf_utils.mod.f90:$(SRCDIR)/crotwf_utils.mod.F90
crotwf_utils.mod.o:crotwf_utils.mod.f90 error_handling.mod \
                kinds.mod mp_interface.mod ovlap_utils.mod \
                parac.mod rotate_utils.mod spin.mod system.mod \
                utils.mod zeroing_utils.mod

csize_utils.mod.f90:$(SRCDIR)/csize_utils.mod.F90
csize_utils.mod.o:csize_utils.mod.f90 dotp_utils.mod elct.mod \
                error_handling.mod kinds.mod kpts.mod mp_interface.mod \
                parac.mod spin.mod system.mod utils.mod

csmat_utils.mod.f90:$(SRCDIR)/csmat_utils.mod.F90
csmat_utils.mod.o:csmat_utils.mod.f90 cvan.mod error_handling.mod \
                ions.mod jrotation_utils.mod kinds.mod mp_interface.mod \
                nlps.mod nort.mod ovlap_utils.mod parac.mod \
                pslo.mod spin.mod system.mod timer.mod

cublas_interfaces.mod.f90:$(SRCDIR)/cublas_interfaces.mod.F90
cublas_interfaces.mod.o:cublas_interfaces.mod.f90 cuda_interfaces.mod

cublas_types.mod.f90:$(SRCDIR)/cublas_types.mod.F90
cublas_types.mod.o:cublas_types.mod.f90 cublas_interfaces.mod \
                cuda_types.mod

cublas_utils.mod.f90:$(SRCDIR)/cublas_utils.mod.F90
cublas_utils.mod.o:cublas_utils.mod.f90 cublas_interfaces.mod \
                cublas_types.mod cuda_types.mod cuda_utils.mod \
                error_handling.mod kinds.mod sizeof_kinds.mod \
                string_utils.mod

cuda_interfaces.mod.f90:$(SRCDIR)/cuda_interfaces.mod.F90
cuda_interfaces.mod.o:cuda_interfaces.mod.f90

cuda_types.mod.f90:$(SRCDIR)/cuda_types.mod.F90
cuda_types.mod.o:cuda_types.mod.f90 cuda_interfaces.mod kinds.mod

cuda_utils.mod.f90:$(SRCDIR)/cuda_utils.mod.F90
cuda_utils.mod.o:cuda_utils.mod.f90 cuda_interfaces.mod cuda_types.mod \
                error_handling.mod kinds.mod sizeof_kinds.mod \
                string_utils.mod

cufft_interfaces.mod.f90:$(SRCDIR)/cufft_interfaces.mod.F90
cufft_interfaces.mod.o:cufft_interfaces.mod.f90 cuda_interfaces.mod \
                kinds.mod

cufft_types.mod.f90:$(SRCDIR)/cufft_types.mod.F90
cufft_types.mod.o:cufft_types.mod.f90 cuda_types.mod cufft_interfaces.mod

cufft_utils.mod.f90:$(SRCDIR)/cufft_utils.mod.F90
cufft_utils.mod.o:cufft_utils.mod.f90 cuda_types.mod cuda_utils.mod \
                cufft_interfaces.mod cufft_types.mod error_handling.mod \
                string_utils.mod kinds.mod

cusolver_interfaces.mod.f90:$(SRCDIR)/cusolver_interfaces.mod.F90
cusolver_interfaces.mod.o:cusolver_interfaces.mod.f90 cublas_interfaces.mod \
                cuda_interfaces.mod

cusolver_types.mod.f90:$(SRCDIR)/cusolver_types.mod.F90
cusolver_types.mod.o:cusolver_types.mod.f90 cuda_types.mod \
                cusolver_interfaces.mod

cusolver_utils.mod.f90:$(SRCDIR)/cusolver_utils.mod.F90
cusolver_utils.mod.o:cusolver_utils.mod.f90 cublas_interfaces.mod \
                cuda_types.mod cuda_utils.mod cusolver_interfaces.mod \
                cusolver_types.mod error_handling.mod kinds.mod \
                sizeof_kinds.mod string_utils.mod

cuuser_interfaces.mod.f90:$(SRCDIR)/cuuser_interfaces.mod.F90
cuuser_interfaces.mod.o:cuuser_interfaces.mod.f90 cuda_interfaces.mod \
                kinds.mod

cuuser_utils.mod.f90:$(SRCDIR)/cuuser_utils.mod.F90
cuuser_utils.mod.o:cuuser_utils.mod.f90 cuda_types.mod cuda_utils.mod \
                cuuser_interfaces.mod error_handling.mod kinds.mod

cvan.mod.f90:   $(SRCDIR)/cvan.mod.F90
cvan.mod.o:     cvan.mod.f90 kinds.mod

davidson_utils.mod.f90:$(SRCDIR)/davidson_utils.mod.F90
davidson_utils.mod.o:davidson_utils.mod.f90 cvan.mod dotp_utils.mod \
                error_handling.mod gsortho_utils.mod hfx_drivers.mod \
                hpsi_utils.mod kinds.mod machine.mod mp_interface.mod \
                parac.mod pslo.mod rnlsm_utils.mod soft.mod \
                sort_utils.mod spsi_utils.mod system.mod testex_utils.mod \
                timer.mod utils.mod vgsortho_utils.mod zeroing_utils.mod

dcacp_utils.mod.f90:$(SRCDIR)/dcacp_utils.mod.F90
dcacp_utils.mod.o:dcacp_utils.mod.f90 atom.mod cp_xc_utils.mod \
                dpot.mod error_handling.mod func.mod ions.mod \
                kinds.mod mp_interface.mod nlps.mod parac.mod \
                pslo.mod recpnew_utils.mod sgpp.mod system.mod \
                timer.mod

dd_functionals_utils.mod.f90:$(SRCDIR)/dd_functionals_utils.mod.F90
dd_functionals_utils.mod.o:dd_functionals_utils.mod.f90 kinds.mod

ddip.mod.f90:   $(SRCDIR)/ddip.mod.F90
ddip.mod.o:     ddip.mod.f90 kinds.mod

ddipo_utils.mod.f90:$(SRCDIR)/ddipo_utils.mod.F90
ddipo_utils.mod.o:ddipo_utils.mod.f90 atwf.mod ddip.mod dipomod.mod \
                dotp_utils.mod elct.mod error_handling.mod \
                forcedr_utils.mod gvec.mod ions.mod jrotation_utils.mod \
                kinds.mod kpts.mod latgen_utils.mod metr.mod \
                mp_interface.mod numpw_utils.mod opeigr_utils.mod \
                parac.mod prcp.mod rggen_utils.mod rotate_utils.mod \
                sd_wannier_utils.mod setbasis_utils.mod sort_utils.mod \
                sphe.mod spin.mod system.mod timer.mod utils.mod \
                wann.mod wannier_center_utils.mod zeroing_utils.mod

dd_xc_ana_utils.mod.f90:$(SRCDIR)/dd_xc_ana_utils.mod.F90
dd_xc_ana_utils.mod.o:dd_xc_ana_utils.mod.f90 cnst.mod cppt.mod \
                cp_xc_utils.mod cp_dxc_driver.mod dd_functionals_utils.mod \
                error_handling.mod fft_maxfft.mod fftmain_utils.mod \
                fftnew_utils.mod func.mod graden_utils.mod \
                kinds.mod lr_xcpot_utils.mod mp_interface.mod \
                parac.mod spin.mod system.mod timer.mod utils.mod \
                zeroing_utils.mod

dd_xc_utils.mod.f90:$(SRCDIR)/dd_xc_utils.mod.F90
dd_xc_utils.mod.o:dd_xc_utils.mod.f90 cppt.mod dd_xc_ana_utils.mod \
                error_handling.mod fftmain_utils.mod gcener_utils.mod \
                graden_utils.mod kinds.mod linres.mod lr_xcpot_utils.mod \
                nlcc.mod parac.mod spin.mod system.mod timer.mod \
                utils.mod xcener_utils.mod zeroing_utils.mod

debfor_utils.mod.f90:$(SRCDIR)/debfor_utils.mod.F90
debfor_utils.mod.o:debfor_utils.mod.f90 adat.mod andp.mod calc_alm_utils.mod \
                coor.mod copot_utils.mod cotr.mod detdof_utils.mod \
                dynit_utils.mod ehpsi_utils.mod elct.mod ener.mod \
                error_handling.mod fint.mod forcedr_driver.mod \
                forcedr_utils.mod forcep_utils.mod initrun_driver.mod \
                initrun_utils.mod ions.mod kinds.mod kpts.mod \
                linres.mod lr_tddft_utils.mod machine.mod nlcc.mod \
                norm.mod ortho_utils.mod parac.mod pbc_utils.mod \
                phfac_utils.mod poin.mod pslo.mod rhoofr_utils.mod \
                rnlsm_utils.mod ropt.mod rrane_utils.mod setirec_utils.mod \
                soft.mod spin.mod store_types.mod symm.mod \
                symtrz_utils.mod system.mod testex_utils.mod \
                updrho_utils.mod updwf_utils.mod wrener_utils.mod \
                wrgeo_utils.mod wv30_utils.mod zeroing_utils.mod

density_functionals_utils.mod.f90:$(SRCDIR)/density_functionals_utils.mod.F90
density_functionals_utils.mod.o:density_functionals_utils.mod.f90 \
                error_handling.mod kinds.mod

density_utils.mod.f90:$(SRCDIR)/density_utils.mod.F90
density_utils.mod.o:density_utils.mod.f90 kinds.mod

densrd_utils.mod.f90:$(SRCDIR)/densrd_utils.mod.F90
densrd_utils.mod.o:densrd_utils.mod.f90 error_handling.mod \
                fileopen_utils.mod fileopenmod.mod gvec.mod \
                ions.mod kinds.mod mp_interface.mod parac.mod \
                response_pmod.mod system.mod

densto_utils.mod.f90:$(SRCDIR)/densto_utils.mod.F90
densto_utils.mod.o:densto_utils.mod.f90 fileopen_utils.mod \
                kinds.mod error_handling.mod timer.mod mp_interface.mod \
                parac.mod fileopen_utils.mod fileopen_utils.mod \
                system.mod kinds.mod error_handling.mod timer.mod \
                mp_interface.mod system.mod parac.mod cell.mod \
                gvec.mod ions.mod fileopenmod.mod fileopen_utils.mod \
                fileopen_utils.mod response_pmod.mod mm_dimmod.mod

deort_utils.mod.f90:$(SRCDIR)/deort_utils.mod.F90
deort_utils.mod.o:deort_utils.mod.f90 dotp_utils.mod error_handling.mod \
                kinds.mod mp_interface.mod parac.mod timer.mod \
                utils.mod

detdof_utils.mod.f90:$(SRCDIR)/detdof_utils.mod.F90
detdof_utils.mod.o:detdof_utils.mod.f90 adat.mod cnstfc_utils.mod \
                cnstpr_utils.mod cotr.mod dum2_utils.mod error_handling.mod \
                ions.mod isos.mod kinds.mod mm_dimmod.mod mm_input.mod \
                nose.mod parac.mod pimd.mod rmas.mod system.mod \
                tpar.mod zeroing_utils.mod

detsp_utils.mod.f90:$(SRCDIR)/detsp_utils.mod.F90
detsp_utils.mod.o:detsp_utils.mod.f90 array_utils.mod atom.mod \
                clas.mod coor.mod cotr.mod error_handling.mod \
                fileopen_utils.mod fileopenmod.mod if_parallel.mod \
                inscan_utils.mod ions.mod kinds.mod mm_input.mod \
                mm_parallel.mod mp_interface.mod nlcc.mod nlps.mod \
                parac.mod prmem_utils.mod readsr_utils.mod \
                sgpp.mod system.mod timer.mod zeroing_utils.mod

dftin_utils.mod.f90:$(SRCDIR)/dftin_utils.mod.F90
dftin_utils.mod.o:dftin_utils.mod.f90 cp_gga_correlation_utils.mod \
                cp_gga_exchange_utils.mod cp_xc_utils.mod hubbardu.mod \
                ener.mod error_handling.mod func.mod hfxmod.mod \
                initclust_utils.mod inscan_utils.mod kinds.mod \
                linres.mod mp_interface.mod parac.mod pw_hfx_input_cnst.mod \
                readsr_utils.mod store_types.mod system.mod \
                tbxc.mod vdwcmod.mod wann.mod zeroing_utils.mod

dginit_utils.mod.f90:$(SRCDIR)/dginit_utils.mod.F90
dginit_utils.mod.o:dginit_utils.mod.f90 dg.mod fftnew_utils.mod \
                timer.mod

dg.mod.f90:     $(SRCDIR)/dg.mod.F90
dg.mod.o:       dg.mod.f90 kinds.mod

difrho_utils.mod.f90:$(SRCDIR)/difrho_utils.mod.F90
difrho_utils.mod.o:difrho_utils.mod.f90 cnst.mod cppt.mod elct.mod \
                error_handling.mod fft.mod fft_maxfft.mod fftmain_utils.mod \
                geq0mod.mod kinds.mod mp_interface.mod parac.mod \
                pslo.mod reshaper.mod rhov_utils.mod spin.mod \
                system.mod timer.mod zeroing_utils.mod

dipo.mod.f90:   $(SRCDIR)/dipo.mod.F90
dipo.mod.o:     dipo.mod.f90 kinds.mod

dipo_utils.mod.f90:$(SRCDIR)/dipo_utils.mod.F90
dipo_utils.mod.o:dipo_utils.mod.f90 bessm_utils.mod cnst.mod \
                cppt.mod dipomod.mod error_handling.mod fft_maxfft.mod \
                fftmain_utils.mod fftnew_utils.mod geq0mod.mod \
                ions.mod kinds.mod mp_interface.mod parac.mod \
                system.mod timer.mod

disortho_utils.mod.f90:$(SRCDIR)/disortho_utils.mod.F90
disortho_utils.mod.o:disortho_utils.mod.f90 cp_cuda_types.mod \
                cp_cuortho_types.mod cp_grp_utils.mod cublas_types.mod \
                cublas_utils.mod cuda_types.mod cuda_utils.mod \
                cusolver_types.mod cuuser_utils.mod error_handling.mod \
                gsortho_utils.mod jrotation_utils.mod kinds.mod \
                linalg_utils.mod mp_interface.mod ovlap_utils.mod \
                parac.mod part_1d.mod string_utils.mod system.mod \
                thread_view_types.mod thread_view_utils.mod \
                timer.mod zeroing_utils.mod nvtx_utils.mod

dispp_utils.mod.f90:$(SRCDIR)/dispp_utils.mod.F90
dispp_utils.mod.o:dispp_utils.mod.f90 ions.mod kinds.mod system.mod

dist_friesner_utils.mod.f90:$(SRCDIR)/dist_friesner_utils.mod.F90
dist_friesner_utils.mod.o:dist_friesner_utils.mod.f90 adjmu_utils.mod \
                dotp_utils.mod error_handling.mod fint.mod \
                friesner_utils.mod gsortho_utils.mod hfx_drivers.mod \
                hpsi_utils.mod jrotation_utils.mod kinds.mod \
                kpts.mod machine.mod mp_interface.mod parac.mod \
                prng_utils.mod randtowf_utils.mod sort_utils.mod \
                summat_utils.mod system.mod timer.mod zeroing_utils.mod

dist_prowfn_utils.mod.f90:$(SRCDIR)/dist_prowfn_utils.mod.F90
dist_prowfn_utils.mod.o:dist_prowfn_utils.mod.f90 adat.mod \
                atwf.mod augchg_utils.mod cmaos_utils.mod dotp_utils.mod \
                elct.mod error_handling.mod fileopen_utils.mod \
                fileopenmod.mod forcep_utils.mod ions.mod jrotation_utils.mod \
                kinds.mod linalg_utils.mod mp_interface.mod \
                parac.mod prden.mod prmem_utils.mod prng_utils.mod \
                prop.mod prowfn_utils.mod pslo.mod rnlsm_utils.mod \
                setbasis_utils.mod sfac.mod spin.mod system.mod \
                timer.mod zeroing_utils.mod

d_mat_p_utils.mod.f90:$(SRCDIR)/d_mat_p_utils.mod.F90
d_mat_p_utils.mod.o:d_mat_p_utils.mod.f90 cnst.mod cppt.mod \
                error_handling.mod geq0mod.mod ions.mod kinds.mod \
                nlps.mod parac.mod pbc_utils.mod pslo.mod ragg.mod \
                sfac.mod sgpp.mod special_functions.mod system.mod \
                timer.mod kinds.mod error_handling.mod timer.mod \
                system.mod parac.mod cppt.mod geq0mod.mod sfac.mod \
                kinds.mod error_handling.mod timer.mod system.mod \
                parac.mod cppt.mod geq0mod.mod

dnlpdk_p_utils.mod.f90:$(SRCDIR)/dnlpdk_p_utils.mod.F90
dnlpdk_p_utils.mod.o:dnlpdk_p_utils.mod.f90 cppt.mod dpot.mod \
                error_handling.mod fint.mod fitpack_utils.mod \
                ions.mod kinds.mod nlcc.mod nlps.mod parac.mod \
                prmem_utils.mod pslo.mod qspl.mod response_pmod.mod \
                sgpp.mod system.mod timer.mod ylmr2_utils.mod \
                zeroing_utils.mod

domdr_utils.mod.f90:$(SRCDIR)/domdr_utils.mod.F90
domdr_utils.mod.o:domdr_utils.mod.f90 cppt.mod cvan.mod hubbardu.mod \
                elct.mod error_handling.mod fnlalloc_utils.mod \
                geq0mod.mod ions.mod kinds.mod kpts.mod mp_interface.mod \
                nlps.mod parac.mod pslo.mod rnlsm_utils.mod \
                sfac.mod spin.mod spsi_utils.mod system.mod \
                timer.mod zeroing_utils.mod

do_perturbation_p_utils.mod.f90:$(SRCDIR)/do_perturbation_p_utils.mod.F90
do_perturbation_p_utils.mod.o:do_perturbation_p_utils.mod.f90 \
                coor.mod cotr.mod detdof_utils.mod eicalc_utils.mod \
                eigensystem_p_utils.mod epr_p_utils.mod error_handling.mod \
                forcedr_utils.mod forcep_utils.mod forces_p_utils.mod \
                fukui_p_utils.mod hardness_p_utils.mod interaction_manno_p_utils.mod \
                interaction_p_utils.mod isos.mod kinds.mod \
                kpts.mod lanc_phon_p_utils.mod machine.mod \
                mm_input.mod mp_interface.mod nmr_p_utils.mod \
                nuclear_p_utils.mod opeigr_utils.mod parac.mod \
                pert_kpoint_p_utils.mod perturbation_p_utils.mod \
                phfac_utils.mod phonons_p_utils.mod prmem_utils.mod \
                pslo.mod raman_p_utils.mod response_pmod.mod \
                rhoofr_p_utils.mod rhoofr_utils.mod rhopri_utils.mod \
                rnlsm_utils.mod ropt.mod rscpot_utils.mod spin.mod \
                store_types.mod symm.mod system.mod updwf_p_utils.mod \
                voa_p_utils.mod vofrho_utils.mod zeroing_utils.mod

dotp_utils.mod.f90:$(SRCDIR)/dotp_utils.mod.F90
dotp_utils.mod.o:dotp_utils.mod.f90 error_handling.mod geq0mod.mod \
                kinds.mod timer.mod kinds.mod error_handling.mod \
                timer.mod geq0mod.mod

dpot.mod.f90:   $(SRCDIR)/dpot.mod.F90
dpot.mod.o:     dpot.mod.f90 system.mod

dqgalloc_utils.mod.f90:$(SRCDIR)/dqgalloc_utils.mod.F90
dqgalloc_utils.mod.o:dqgalloc_utils.mod.f90 error_handling.mod \
                fft_maxfft.mod pslo.mod spin.mod str2.mod system.mod

dqvan2_utils.mod.f90:$(SRCDIR)/dqvan2_utils.mod.F90
dqvan2_utils.mod.o:dqvan2_utils.mod.f90 zeroing_utils.mod kinds.mod \
                error_handling.mod timer.mod system.mod parac.mod \
                nlps.mod cvan.mod geq0mod.mod qspl.mod cppt.mod \
                aavan.mod dylmr_utils.mod ylmr2_utils.mod fitpack_utils.mod \
                zeroing_utils.mod

drhov_utils.mod.f90:$(SRCDIR)/drhov_utils.mod.F90
drhov_utils.mod.o:drhov_utils.mod.f90 zeroing_utils.mod kinds.mod \
                error_handling.mod timer.mod mp_interface.mod \
                prmem_utils.mod reshaper.mod system.mod parac.mod \
                ions.mod pslo.mod nlps.mod elct.mod cppt.mod \
                sfac.mod spin.mod str2.mod zeroing_utils.mod

dum2_utils.mod.f90:$(SRCDIR)/dum2_utils.mod.F90
dum2_utils.mod.o:dum2_utils.mod.f90 cotr.mod error_handling.mod \
                ions.mod kinds.mod mm_dimmod.mod parac.mod \
                rmas.mod system.mod

dylmr_utils.mod.f90:$(SRCDIR)/dylmr_utils.mod.F90
dylmr_utils.mod.o:dylmr_utils.mod.f90 cnst.mod error_handling.mod \
                geq0mod.mod kinds.mod parac.mod str2.mod strs.mod \
                system.mod zeroing_utils.mod

dynit_utils.mod.f90:$(SRCDIR)/dynit_utils.mod.F90
dynit_utils.mod.o:dynit_utils.mod.f90 clas.mod cnst.mod ions.mod \
                kinds.mod rmas.mod system.mod tpar.mod

eam.mod.f90:    $(SRCDIR)/eam.mod.F90
eam.mod.o:      eam.mod.f90 kinds.mod system.mod

eam_pot_utils.mod.f90:$(SRCDIR)/eam_pot_utils.mod.F90
eam_pot_utils.mod.o:eam_pot_utils.mod.f90 adat.mod cnst.mod \
                dpot.mod eam.mod error_handling.mod inscan_utils.mod \
                ions.mod kinds.mod metr.mod mp_interface.mod \
                parac.mod pbc_utils.mod system.mod timer.mod \
                zeroing_utils.mod

eextern_utils.mod.f90:$(SRCDIR)/eextern_utils.mod.F90
eextern_utils.mod.o:eextern_utils.mod.f90 atomc_utils.mod cnst.mod \
                efld.mod epot_types.mod ions.mod kinds.mod \
                metr.mod mm_input.mod parac.mod pbc_utils.mod \
                ragg.mod system.mod timer.mod

efield_utils.mod.f90:$(SRCDIR)/efield_utils.mod.F90
efield_utils.mod.o:efield_utils.mod.f90 ddip.mod ddipo_utils.mod \
                ener.mod error_handling.mod geq0mod.mod gvec.mod \
                ions.mod kinds.mod opeigr_p_utils.mod parac.mod \
                system.mod timer.mod utils.mod zeroing_utils.mod

efld.mod.f90:   $(SRCDIR)/efld.mod.F90
efld.mod.o:     efld.mod.f90 kinds.mod

egointer_utils.mod.f90:$(SRCDIR)/egointer_utils.mod.F90
egointer_utils.mod.o:egointer_utils.mod.f90 atomc_utils.mod \
                atwf.mod cnst.mod coor.mod copot_utils.mod \
                cppt.mod difrho_utils.mod dipo_utils.mod dipomod.mod \
                dynit_utils.mod efld.mod eicalc_utils.mod elct.mod \
                elstpo_utils.mod ener.mod epot_types.mod error_handling.mod \
                espchg_utils.mod exdipo_utils.mod fft_maxfft.mod \
                fftmain_utils.mod fftnew_utils.mod fileopenmod.mod \
                finalp_utils.mod forcedr_driver.mod forcedr_utils.mod \
                forcep_utils.mod geq0mod.mod hip_utils.mod \
                inscan_utils.mod ions.mod isos.mod kinds.mod \
                kpts.mod lodipo_utils.mod lodp.mod machine.mod \
                mp_interface.mod mulliken_utils.mod nlcc.mod \
                norm.mod parac.mod phfac_utils.mod prop.mod \
                propin_utils.mod proppt_utils.mod pslo.mod \
                purge_utils.mod readsr_utils.mod rhoofr_c_utils.mod \
                rhoofr_utils.mod rhopri_utils.mod rnlsm_utils.mod \
                ropt.mod rrane_utils.mod setirec_utils.mod \
                special_functions.mod spin.mod store_types.mod \
                system.mod testex_utils.mod timer.mod updwf_utils.mod \
                wrener_utils.mod wrgeo_utils.mod wv30_utils.mod \
                zeroing_utils.mod

ehpsi_utils.mod.f90:$(SRCDIR)/ehpsi_utils.mod.F90
ehpsi_utils.mod.o:ehpsi_utils.mod.f90 cppt.mod error_handling.mod \
                fft.mod fftmain_utils.mod fint.mod geq0mod.mod \
                ions.mod kinds.mod kpclean_utils.mod kpnt.mod \
                kpts.mod mp_interface.mod nlps.mod parac.mod \
                reshaper.mod sgpp.mod spin.mod system.mod timer.mod \
                zeroing_utils.mod zeroing_utils.mod cnst.mod \
                cppt.mod ehpsi_utils.mod error_handling.mod \
                fft.mod fft_maxfft.mod fftmain_utils.mod fint.mod \
                geq0mod.mod kinds.mod mp_interface.mod nlps.mod \
                parac.mod reshaper.mod system.mod timer.mod

ehrenfest_utils.mod.f90:$(SRCDIR)/ehrenfest_utils.mod.F90
ehrenfest_utils.mod.o:ehrenfest_utils.mod.f90 elct.mod error_handling.mod \
                kinds.mod mp_interface.mod parac.mod spin.mod \
                system.mod td_cayley_utils.mod td_input.mod \
                timer.mod zeroing_utils.mod

eicalc_utils.mod.f90:$(SRCDIR)/eicalc_utils.mod.F90
eicalc_utils.mod.o:eicalc_utils.mod.f90 cppt.mod ions.mod kinds.mod \
                nvtx_utils.mod sfac.mod system.mod timer.mod \
                zeroing_utils.mod

eigensystem_p_utils.mod.f90:$(SRCDIR)/eigensystem_p_utils.mod.F90
eigensystem_p_utils.mod.o:eigensystem_p_utils.mod.f90 coor.mod \
                cppt.mod dotp_utils.mod error_handling.mod \
                fft_maxfft.mod fftmain_utils.mod kinds.mod \
                mp_interface.mod parac.mod readsr_utils.mod \
                response_pmod.mod rhoofr_p_utils.mod rhoofr_utils.mod \
                rwfopt_p_utils.mod spin.mod system.mod utils.mod

eind_ii_utils.mod.f90:$(SRCDIR)/eind_ii_utils.mod.F90
eind_ii_utils.mod.o:eind_ii_utils.mod.f90 cnst.mod ions.mod \
                kinds.mod metr.mod parac.mod pbc_utils.mod \
                ragg.mod special_functions.mod system.mod timer.mod

eind_loc_utils.mod.f90:$(SRCDIR)/eind_loc_utils.mod.F90
eind_loc_utils.mod.o:eind_loc_utils.mod.f90 cppt.mod fftmain_utils.mod \
                geq0mod.mod kinds.mod parac.mod sfac.mod system.mod

eind_nl_utils.mod.f90:$(SRCDIR)/eind_nl_utils.mod.F90
eind_nl_utils.mod.o:eind_nl_utils.mod.f90 elct.mod error_handling.mod \
                kinds.mod nlps.mod parac.mod pslo.mod sfac.mod \
                sgpp.mod system.mod

ekinpp_utils.mod.f90:$(SRCDIR)/ekinpp_utils.mod.F90
ekinpp_utils.mod.o:ekinpp_utils.mod.f90 ions.mod kinds.mod \
                pimd.mod rmas.mod timer.mod

elct2.mod.f90:  $(SRCDIR)/elct2.mod.F90
elct2.mod.o:    elct2.mod.f90

elct.mod.f90:   $(SRCDIR)/elct.mod.F90
elct.mod.o:     elct.mod.f90 kinds.mod

elf_utils.mod.f90:$(SRCDIR)/elf_utils.mod.F90
elf_utils.mod.o:elf_utils.mod.f90 bsym.mod cnst.mod cppt.mod \
                elct.mod error_handling.mod fft_maxfft.mod \
                fftmain_utils.mod fftnew_utils.mod kinds.mod \
                kpts.mod meta_multiple_walkers_utils.mod mp_interface.mod \
                mw.mod parac.mod phfac_utils.mod pimd.mod prden.mod \
                pslo.mod readsr_utils.mod rhoofr_c_utils.mod \
                rhoofr_utils.mod rnlsm_utils.mod spin.mod system.mod \
                zeroing_utils.mod

elstpo_utils.mod.f90:$(SRCDIR)/elstpo_utils.mod.F90
elstpo_utils.mod.o:elstpo_utils.mod.f90 cppt.mod error_handling.mod \
                fftmain_utils.mod fftnew_utils.mod geq0mod.mod \
                hip_utils.mod isos.mod kinds.mod parac.mod \
                system.mod timer.mod zeroing_utils.mod

empf.mod.f90:   $(SRCDIR)/empf.mod.F90
empf.mod.o:     empf.mod.f90 kinds.mod

empfor_utils.mod.f90:$(SRCDIR)/empfor_utils.mod.F90
empfor_utils.mod.o:empfor_utils.mod.f90 adat.mod cnst.mod cotr.mod \
                empf.mod error_handling.mod fstart_utils.mod \
                ions.mod kinds.mod parac.mod timer.mod zeroing_utils.mod

enbandpri_utils.mod.f90:$(SRCDIR)/enbandpri_utils.mod.F90
enbandpri_utils.mod.o:enbandpri_utils.mod.f90 cnst.mod elct.mod \
                fileopen_utils.mod fileopenmod.mod kinds.mod \
                kpnt.mod parac.mod system.mod timer.mod

ener.mod.f90:   $(SRCDIR)/ener.mod.F90
ener.mod.o:     ener.mod.f90 kinds.mod system.mod

enosmove_utils.mod.f90:$(SRCDIR)/enosmove_utils.mod.F90
enosmove_utils.mod.o:enosmove_utils.mod.f90 kinds.mod nose.mod \
                system.mod

envir_utils.mod.f90:$(SRCDIR)/envir_utils.mod.F90
envir_utils.mod.o:envir_utils.mod.f90 envj.mod machine.mod

envj.mod.f90:   $(SRCDIR)/envj.mod.F90
envj.mod.o:     envj.mod.f90 kinds.mod

epot_types.mod.f90:$(SRCDIR)/epot_types.mod.F90
epot_types.mod.o:epot_types.mod.f90 kinds.mod

epr_current_p_utils.mod.f90:$(SRCDIR)/epr_current_p_utils.mod.F90
epr_current_p_utils.mod.o:epr_current_p_utils.mod.f90 coor.mod \
                error_handling.mod forcep_utils.mod geq0mod.mod \
                ions.mod kinds.mod nmr_position_p_utils.mod \
                nmr_util_p_utils.mod parac.mod prmem_utils.mod \
                response_pmod.mod spin.mod system.mod timer.mod \
                zeroing_utils.mod

epr_dv0_utils.mod.f90:$(SRCDIR)/epr_dv0_utils.mod.F90
epr_dv0_utils.mod.o:epr_dv0_utils.mod.f90 cnst.mod cppt.mod \
                fft_maxfft.mod fftmain_utils.mod fftnew_utils.mod \
                kinds.mod parac.mod system.mod timer.mod utils.mod \
                zeroing_utils.mod

epr_efg_utils.mod.f90:$(SRCDIR)/epr_efg_utils.mod.F90
epr_efg_utils.mod.o:epr_efg_utils.mod.f90 cnst.mod cppt.mod \
                eicalc_utils.mod error_handling.mod fileopen_utils.mod \
                fileopenmod.mod geq0mod.mod ions.mod kinds.mod \
                parac.mod prop.mod sfac.mod spin.mod system.mod \
                timer.mod zeroing_utils.mod kinds.mod error_handling.mod \
                timer.mod mp_interface.mod system.mod parac.mod \
                ions.mod spin.mod epr_efg_utils.mod eicalc_utils.mod \
                zeroing_utils.mod

epr_hyp_utils.mod.f90:$(SRCDIR)/epr_hyp_utils.mod.F90
epr_hyp_utils.mod.o:epr_hyp_utils.mod.f90 adat.mod atwf.mod \
                cnst.mod cppt.mod error_handling.mod fftmain_utils.mod \
                fitpack_utils.mod geq0mod.mod ions.mod kinds.mod \
                mp_interface.mod parac.mod pslo.mod sfac.mod \
                system.mod zeroing_utils.mod

epr_p_utils.mod.f90:$(SRCDIR)/epr_p_utils.mod.F90
epr_p_utils.mod.o:epr_p_utils.mod.f90 cnst.mod coor.mod cppt.mod \
                ddip.mod ddipo_utils.mod eicalc_utils.mod elct.mod \
                epr_current_p_utils.mod epr_dv0_utils.mod epr_hyp_utils.mod \
                epr_util_p_utils.mod error_handling.mod fft_maxfft.mod \
                fftmain_utils.mod fftnew_utils.mod fitpack_utils.mod \
                geq0mod.mod ions.mod isos.mod kinds.mod localize_utils.mod \
                machine.mod mp_interface.mod nmr_full_p_utils.mod \
                nmr_position_p_utils.mod nmr_util_p_utils.mod \
                parac.mod perturbation_p_utils.mod phfac_utils.mod \
                prmem_utils.mod prop.mod qspl.mod ragg.mod \
                response_pmod.mod restart_p_utils.mod rwfopt_p_utils.mod \
                sfac.mod soft.mod special_functions.mod spin.mod \
                store_types.mod system.mod timer.mod utils.mod \
                wann.mod xcener_utils.mod zeroing_utils.mod

epr_util_p_utils.mod.f90:$(SRCDIR)/epr_util_p_utils.mod.F90
epr_util_p_utils.mod.o:epr_util_p_utils.mod.f90 error_handling.mod \
                kinds.mod nmr_position_p_utils.mod nmr_util_p_utils.mod \
                parac.mod perturbation_p_utils.mod response_pmod.mod \
                rwfopt_p_utils.mod soft.mod system.mod timer.mod \
                zeroing_utils.mod

error_handling.mod.f90:$(SRCDIR)/error_handling.mod.F90
error_handling.mod.o:error_handling.mod.f90 parac.mod

espchg_utils.mod.f90:$(SRCDIR)/espchg_utils.mod.F90
espchg_utils.mod.o:espchg_utils.mod.f90 adat.mod cnst.mod cppt.mod \
                elct.mod error_handling.mod fftmain_utils.mod \
                fftnew_utils.mod geq0mod.mod hip_utils.mod \
                ions.mod isos.mod kinds.mod mp_interface.mod \
                parac.mod pbc_utils.mod sfac.mod system.mod \
                timer.mod zeroing_utils.mod

evirial_utils.mod.f90:$(SRCDIR)/evirial_utils.mod.F90
evirial_utils.mod.o:evirial_utils.mod.f90 cnst.mod ions.mod \
                kinds.mod nose.mod pbc_utils.mod pimd.mod system.mod

exdipo_utils.mod.f90:$(SRCDIR)/exdipo_utils.mod.F90
exdipo_utils.mod.o:exdipo_utils.mod.f90 dipo_utils.mod dipomod.mod \
                kinds.mod lodipo_utils.mod lodp.mod prop.mod \
                system.mod timer.mod

exterp_utils.mod.f90:$(SRCDIR)/exterp_utils.mod.F90
exterp_utils.mod.o:exterp_utils.mod.f90 efld.mod error_handling.mod \
                extpotmod.mod isos.mod kinds.mod mp_interface.mod \
                parac.mod prmem_utils.mod system.mod timer.mod \
                zeroing_utils.mod

extpot.mod.f90: $(SRCDIR)/extpot.mod.F90
extpot.mod.o:   extpot.mod.f90 kinds.mod

extrap_utils.mod.f90:$(SRCDIR)/extrap_utils.mod.F90
extrap_utils.mod.o:extrap_utils.mod.f90 kinds.mod

fcas.mod.f90:   $(SRCDIR)/fcas.mod.F90
fcas.mod.o:     fcas.mod.f90 kinds.mod

ffsum_utils.mod.f90:$(SRCDIR)/ffsum_utils.mod.F90
ffsum_utils.mod.o:ffsum_utils.mod.f90 cppt.mod fitpack_utils.mod \
                geq0mod.mod ions.mod kinds.mod qspl.mod sfac.mod \
                system.mod

fftchk_utils.mod.f90:$(SRCDIR)/fftchk_utils.mod.F90
fftchk_utils.mod.o:fftchk_utils.mod.f90 cp_cuda_types.mod error_handling.mod \
                parac.mod

fftcu_methods.mod.f90:$(SRCDIR)/fftcu_methods.mod.F90
fftcu_methods.mod.o:fftcu_methods.mod.f90 cp_cufft_types.mod \
                cp_cufft_utils.mod cublas_types.mod cuda_types.mod \
                cuda_utils.mod cufft_types.mod cuuser_utils.mod \
                fft.mod fft_maxfft.mod fftutil_utils.mod kinds.mod \
                mltfft_utils.mod

fftmain_utils.mod.f90:$(SRCDIR)/fftmain_utils.mod.F90
fftmain_utils.mod.o:fftmain_utils.mod.f90 cp_cuda_types.mod \
                cp_cufft_types.mod cp_cufft_utils.mod cublas_types.mod \
                cuda_types.mod cuda_utils.mod cufft_types.mod \
                error_handling.mod fft.mod fft_maxfft.mod fftcu_methods.mod \
                fftutil_utils.mod kinds.mod mltfft_utils.mod \
                parac.mod system.mod thread_view_types.mod \
                timer.mod

fft_maxfft.mod.f90:$(SRCDIR)/fft_maxfft.mod.F90
fft_maxfft.mod.o:fft_maxfft.mod.f90

fft.mod.f90:    $(SRCDIR)/fft.mod.F90
fft.mod.o:      fft.mod.f90 kinds.mod

fftnew_utils.mod.f90:$(SRCDIR)/fftnew_utils.mod.F90
fftnew_utils.mod.o:fftnew_utils.mod.f90 cell.mod cnst.mod cp_cuda_types.mod \
                cp_cufft_types.mod cppt.mod cuda_types.mod \
                cuda_utils.mod error_handling.mod fft.mod fft_maxfft.mod \
                fftchk_utils.mod kinds.mod loadpa_utils.mod \
                mp_interface.mod parac.mod system.mod utils.mod \
                zeroing_utils.mod

fftprp_utils.mod.f90:$(SRCDIR)/fftprp_utils.mod.F90
fftprp_utils.mod.o:fftprp_utils.mod.f90 cp_cuda_types.mod cp_cufft_types.mod \
                cp_cufft_utils.mod cp_curho_types.mod cp_curho_utils.mod \
                cppt.mod cuda_utils.mod elct.mod error_handling.mod \
                fft.mod fft_maxfft.mod fftnew_utils.mod isos.mod \
                kinds.mod kpts.mod mp_interface.mod parac.mod \
                prmem_utils.mod reshaper.mod rswfmod.mod system.mod \
                timer.mod utils.mod zeroing_utils.mod

fft_utils.mod.f90:$(SRCDIR)/fft_utils.mod.F90
fft_utils.mod.o:fft_utils.mod.f90 error_handling.mod fft.mod

fftutil_utils.mod.f90:$(SRCDIR)/fftutil_utils.mod.F90
fftutil_utils.mod.o:fftutil_utils.mod.f90 fft_maxfft.mod kinds.mod \
                mp_interface.mod parac.mod reshaper.mod system.mod \
                timer.mod utils.mod zeroing_utils.mod

fharm_utils.mod.f90:$(SRCDIR)/fharm_utils.mod.F90
fharm_utils.mod.o:fharm_utils.mod.f90 cnst.mod ions.mod kinds.mod \
                pbc_utils.mod pimd.mod system.mod

fileopen.mod.f90:$(SRCDIR)/fileopen.mod.F90
fileopen.mod.o: fileopen.mod.f90

fileopen_utils.mod.f90:$(SRCDIR)/fileopen_utils.mod.F90
fileopen_utils.mod.o:fileopen_utils.mod.f90 error_handling.mod \
                fileopenmod.mod parac.mod readsr_utils.mod

fillc_utils.mod.f90:$(SRCDIR)/fillc_utils.mod.F90
fillc_utils.mod.o:fillc_utils.mod.f90 cotr.mod ions.mod kinds.mod \
                system.mod

filn.mod.f90:   $(SRCDIR)/filn.mod.F90
filn.mod.o:     filn.mod.f90

finalp_utils.mod.f90:$(SRCDIR)/finalp_utils.mod.F90
finalp_utils.mod.o:finalp_utils.mod.f90 cdft_utils.mod cnst.mod \
                cnstpr_utils.mod dipomod.mod dum2_utils.mod \
                elct.mod ener.mod ions.mod kinds.mod kpts.mod \
                metr.mod norm.mod parac.mod rmas.mod ropt.mod \
                store_types.mod strs.mod struc_utils.mod system.mod \
                timer.mod wrener_utils.mod wrgeo_utils.mod

fint.mod.f90:   $(SRCDIR)/fint.mod.F90
fint.mod.o:     fint.mod.f90 kinds.mod

fitpack_utils.mod.f90:$(SRCDIR)/fitpack_utils.mod.F90
fitpack_utils.mod.o:fitpack_utils.mod.f90 kinds.mod

fixcom_utils.mod.f90:$(SRCDIR)/fixcom_utils.mod.F90
fixcom_utils.mod.o:fixcom_utils.mod.f90 cotr.mod ions.mod kinds.mod

fm_cnst.mod.f90:$(SRCDIR)/fm_cnst.mod.F90
fm_cnst.mod.o:  fm_cnst.mod.f90 cnst.mod kinds.mod

fnlalloc_utils.mod.f90:$(SRCDIR)/fnlalloc_utils.mod.F90
fnlalloc_utils.mod.o:fnlalloc_utils.mod.f90 error_handling.mod \
                ions.mod kinds.mod nlps.mod parac.mod sfac.mod \
                system.mod zeroing_utils.mod

fnonloc_p_utils.mod.f90:$(SRCDIR)/fnonloc_p_utils.mod.F90
fnonloc_p_utils.mod.o:fnonloc_p_utils.mod.f90 cppt.mod error_handling.mod \
                geq0mod.mod kinds.mod kpts.mod mp_interface.mod \
                nlps.mod parac.mod pslo.mod response_pmod.mod \
                sfac.mod sgpp.mod system.mod timer.mod zeroing_utils.mod

fnonloc_utils.mod.f90:$(SRCDIR)/fnonloc_utils.mod.F90
fnonloc_utils.mod.o:fnonloc_utils.mod.f90 cp_grp_utils.mod \
                cppt.mod cvan.mod ener.mod error_handling.mod \
                ions.mod kinds.mod kpnt.mod kpts.mod mp_interface.mod \
                nlps.mod nvtx_utils.mod parac.mod pslo.mod \
                reshaper.mod sfac.mod sgpp.mod spin.mod system.mod \
                timer.mod zeroing_utils.mod

forcedr_driver.mod.f90:$(SRCDIR)/forcedr_driver.mod.F90
forcedr_driver.mod.o:forcedr_driver.mod.f90 bsym.mod cnstfc_utils.mod \
                cotr.mod error_handling.mod forces_driver.mod \
                kinds.mod mm_dim_utils.mod mm_dimmod.mod mp_interface.mod \
                noforce_utils.mod parac.mod pimd.mod reshaper.mod \
                system.mod timer.mod tpar.mod utils.mod zeroing_utils.mod

forcedr_utils.mod.f90:$(SRCDIR)/forcedr_utils.mod.F90
forcedr_utils.mod.o:forcedr_utils.mod.f90 forces_utils.mod \
                k_forces_utils.mod kpts.mod mm_dim_utils.mod \
                mm_dimmod.mod noforce_utils.mod symtrz_utils.mod \
                system.mod timer.mod

forcematch_kfit_utils.mod.f90:$(SRCDIR)/forcematch_kfit_utils.mod.F90
forcematch_kfit_utils.mod.o:forcematch_kfit_utils.mod.f90 error_handling.mod \
                fm_cnst.mod forcematch.mod forcematch_utils.mod \
                kinds.mod mm_dimmod.mod mm_input.mod parac.mod \
                zeroing_utils.mod coordsz.mod $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h $(MODDIR)/Gromos/box.h \
                $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h $(MODDIR)/Gromos/box.h \
                $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h $(MODDIR)/Gromos/box.h \
                $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h $(MODDIR)/Gromos/box.h \
                $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h $(MODDIR)/Gromos/box.h

forcematch.mod.f90:$(SRCDIR)/forcematch.mod.F90
forcematch.mod.o:forcematch.mod.f90 kinds.mod

forcematch_qfit_utils.mod.f90:$(SRCDIR)/forcematch_qfit_utils.mod.F90
forcematch_qfit_utils.mod.o:forcematch_qfit_utils.mod.f90 error_handling.mod \
                fileopen_utils.mod fileopenmod.mod forcematch.mod \
                kinds.mod mm_dimmod.mod parac.mod readsr_utils.mod \
                zeroing_utils.mod coordsz.mod $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h $(MODDIR)/Gromos/box.h

forcematch_utils.mod.f90:$(SRCDIR)/forcematch_utils.mod.F90
forcematch_utils.mod.o:forcematch_utils.mod.f90 error_handling.mod \
                fileopen_utils.mod fileopenmod.mod forcematch.mod \
                kinds.mod machine.mod mm_dimmod.mod parac.mod \
                readsr_utils.mod coordsz.mod $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h $(MODDIR)/Gromos/box.h

forcep_utils.mod.f90:$(SRCDIR)/forcep_utils.mod.F90
forcep_utils.mod.o:forcep_utils.mod.f90 fft_maxfft.mod isos.mod \
                linres.mod spin.mod system.mod

forces_diag_utils.mod.f90:$(SRCDIR)/forces_diag_utils.mod.F90
forces_diag_utils.mod.o:forces_diag_utils.mod.f90 andp.mod \
                andr.mod calc_alm_utils.mod cdft_utils.mod \
                cdftmod.mod cnstfc_utils.mod cotr.mod ehpsi_utils.mod \
                elct.mod elct2.mod ener.mod fint.mod forcedr_driver.mod \
                k_updwf_utils.mod kinds.mod kpts.mod machine.mod \
                mm_extrap.mod mp_interface.mod norm.mod parac.mod \
                pslo.mod rhoofr_c_utils.mod rhoofr_utils.mod \
                rinitwf_driver.mod rnlsm_utils.mod ropt.mod \
                rrandd_utils.mod soft.mod spin.mod store_types.mod \
                system.mod testex_utils.mod timer.mod tpar.mod \
                updrho_utils.mod updwf_utils.mod wrener_utils.mod \
                wv30_utils.mod zeroing_utils.mod

forces_driver.mod.f90:$(SRCDIR)/forces_driver.mod.F90
forces_driver.mod.o:forces_driver.mod.f90 cp_grp_utils.mod \
                csize_utils.mod dotp_utils.mod efield_utils.mod \
                elct.mod ener.mod error_handling.mod fft_maxfft.mod \
                fnonloc_utils.mod func.mod geq0mod.mod gsize_utils.mod \
                hfx_drivers.mod hfxmod.mod hnlmat_utils.mod \
                hubbardu.mod hubbardu_utils.mod jrotation_utils.mod \
                kinds.mod kpts.mod mp_interface.mod nlforce_utils.mod \
                nlps.mod norm.mod nvtx_utils.mod opeigr_utils.mod \
                ovlap_utils.mod parac.mod pslo.mod puttau_utils.mod \
                reshaper.mod rnlfl_utils.mod rnlsm_utils.mod \
                ropt.mod rotate_utils.mod rscpot_utils.mod \
                rswfmod.mod spin.mod summat_utils.mod symtrz_utils.mod \
                system.mod timer.mod utils.mod vpsi_utils.mod \
                zeroing_utils.mod

forces_prop_utils.mod.f90:$(SRCDIR)/forces_prop_utils.mod.F90
forces_prop_utils.mod.o:forces_prop_utils.mod.f90 cnstfc_utils.mod \
                cotr.mod efld.mod ehpsi_utils.mod ehrenfest_utils.mod \
                elct.mod ener.mod error_handling.mod geq0mod.mod \
                k_updwf_utils.mod kinds.mod kpts.mod machine.mod \
                mp_interface.mod nlforce_utils.mod parac.mod \
                ropt.mod summat_utils.mod system.mod td_input.mod \
                timer.mod tpar.mod updrho_utils.mod updwf_utils.mod \
                utils.mod zeroing_utils.mod

forces_p_utils.mod.f90:$(SRCDIR)/forces_p_utils.mod.F90
forces_p_utils.mod.o:forces_p_utils.mod.f90 elct.mod ener.mod \
                h0psi1_p_utils.mod kinds.mod mp_interface.mod \
                norm.mod parac.mod perturbation_p_utils.mod \
                response_pmod.mod rhoofr_p_utils.mod rotate_utils.mod \
                spin.mod system.mod timer.mod v1ofrho_p_utils.mod \
                vpsi_p_utils.mod zeroing_utils.mod

forces_utils.mod.f90:$(SRCDIR)/forces_utils.mod.F90
forces_utils.mod.o:forces_utils.mod.f90 fft_maxfft.mod fnonloc_utils.mod \
                jrotation_utils.mod nlforce_utils.mod nlps.mod \
                opeigr_utils.mod parac.mod pslo.mod rnlsm_utils.mod \
                ropt.mod rscpot_utils.mod spin.mod summat_utils.mod \
                symtrz_utils.mod system.mod

formf_utils.mod.f90:$(SRCDIR)/formf_utils.mod.F90
formf_utils.mod.o:formf_utils.mod.f90 atom.mod cnst.mod dpot.mod \
                error_handling.mod fitpack_utils.mod ions.mod \
                kinds.mod mp_interface.mod parac.mod pslo.mod \
                qspl.mod radin_utils.mod ragg.mod sgpp.mod \
                special_functions.mod system.mod timer.mod \
                utils.mod vdbp.mod zeroing_utils.mod

freqs_utils.mod.f90:$(SRCDIR)/freqs_utils.mod.F90
freqs_utils.mod.o:freqs_utils.mod.f90 cppt.mod cvan.mod elct.mod \
                error_handling.mod harm.mod ions.mod kinds.mod \
                metr.mod mp_interface.mod nlps.mod parac.mod \
                pslo.mod rggen_utils.mod sgpp.mod simulmod.mod \
                system.mod tpar.mod

friesner_c_p_utils.mod.f90:$(SRCDIR)/friesner_c_p_utils.mod.F90
friesner_c_p_utils.mod.o:friesner_c_p_utils.mod.f90 ehpsi_utils.mod \
                error_handling.mod fft_maxfft.mod fint.mod \
                frsblk_c_utils.mod gsortho_utils.mod hpsi_utils.mod \
                kinds.mod parac.mod sort_utils.mod system.mod \
                timer.mod zeroing_utils.mod

friesner_c_utils.mod.f90:$(SRCDIR)/friesner_c_utils.mod.F90
friesner_c_utils.mod.o:friesner_c_utils.mod.f90 adjmu_utils.mod \
                ehpsi_utils.mod error_handling.mod fft_maxfft.mod \
                fint.mod frsblk_c_utils.mod gsortho_utils.mod \
                hpsi_utils.mod jacobi_utils.mod kinds.mod machine.mod \
                mp_interface.mod parac.mod randtowf_utils.mod \
                sort_utils.mod system.mod timer.mod zeroing_utils.mod

friesner_utils.mod.f90:$(SRCDIR)/friesner_utils.mod.F90
friesner_utils.mod.o:friesner_utils.mod.f90 adjmu_utils.mod \
                dotp_utils.mod error_handling.mod fint.mod \
                frsblk_utils.mod gsortho_utils.mod hfx_drivers.mod \
                hpsi_utils.mod jacobi_utils.mod kinds.mod kpts.mod \
                machine.mod mp_interface.mod parac.mod randtowf_utils.mod \
                sort_utils.mod summat_utils.mod system.mod \
                timer.mod zeroing_utils.mod

frsblk_c_utils.mod.f90:$(SRCDIR)/frsblk_c_utils.mod.F90
frsblk_c_utils.mod.o:frsblk_c_utils.mod.f90 adjmu_utils.mod \
                ehpsi_utils.mod error_handling.mod fft_maxfft.mod \
                fint.mod geq0mod.mod gsortho_utils.mod hpsi_utils.mod \
                kinds.mod machine.mod mp_interface.mod parac.mod \
                prng_utils.mod rgs_utils.mod summat_utils.mod \
                system.mod timer.mod zeroing_utils.mod

frsblk_utils.mod.f90:$(SRCDIR)/frsblk_utils.mod.F90
frsblk_utils.mod.o:frsblk_utils.mod.f90 adjmu_utils.mod dotp_utils.mod \
                error_handling.mod fint.mod geq0mod.mod gsortho_utils.mod \
                hpsi_utils.mod kinds.mod kpts.mod machine.mod \
                mp_interface.mod parac.mod prng_utils.mod rgs_utils.mod \
                summat_utils.mod system.mod timer.mod zeroing_utils.mod \
                kinds.mod error_handling.mod timer.mod frsblk_utils.mod \
                hpsi_utils.mod fft_maxfft.mod zeroing_utils.mod \
                kinds.mod error_handling.mod timer.mod frsblk_utils.mod \
                hpsi_utils.mod summat_utils.mod fft_maxfft.mod \
                zeroing_utils.mod kinds.mod error_handling.mod \
                timer.mod rgs_utils.mod

fstart_utils.mod.f90:$(SRCDIR)/fstart_utils.mod.F90
fstart_utils.mod.o:fstart_utils.mod.f90 cnst.mod empf.mod error_handling.mod \
                kinds.mod parac.mod pbc_utils.mod system.mod

fukui_p_utils.mod.f90:$(SRCDIR)/fukui_p_utils.mod.F90
fukui_p_utils.mod.o:fukui_p_utils.mod.f90 cnst.mod coor.mod \
                cppt.mod d_mat_p_utils.mod densrd_utils.mod \
                elct.mod error_handling.mod fft_maxfft.mod \
                fftmain_utils.mod forcep_utils.mod geq0mod.mod \
                ions.mod kinds.mod kpnt.mod mp_interface.mod \
                nlps.mod ovlap_utils.mod parac.mod readsr_utils.mod \
                response_pmod.mod rhoofr_utils.mod rnlfor_utils.mod \
                rnlsm_utils.mod rwfopt_p_utils.mod sfac.mod \
                spin.mod system.mod timer.mod v1ofrho_p_utils.mod \
                vpsi_p_utils.mod wrgeo_utils.mod zeroing_utils.mod

func.mod.f90:   $(SRCDIR)/func.mod.F90
func.mod.o:     func.mod.f90 kinds.mod

functionals_utils.mod.f90:$(SRCDIR)/functionals_utils.mod.F90
functionals_utils.mod.o:functionals_utils.mod.f90 cnst.mod \
                error_handling.mod func.mod kinds.mod special_functions.mod \
                system.mod

fusion_utils.mod.f90:$(SRCDIR)/fusion_utils.mod.F90
fusion_utils.mod.o:fusion_utils.mod.f90 coor.mod elct.mod error_handling.mod \
                filnmod.mod kinds.mod mm_dim_utils.mod mm_dimmod.mod \
                parac.mod readsr_utils.mod rv30_utils.mod setirec_utils.mod \
                shop.mod store_types.mod system.mod wv30_utils.mod \
                zeroing_utils.mod

gcener_utils.mod.f90:$(SRCDIR)/gcener_utils.mod.F90
gcener_utils.mod.o:gcener_utils.mod.f90 cnst.mod cp_xc_driver.mod \
                cppt.mod error_handling.mod fftmain_utils.mod \
                fftnew_utils.mod func.mod functionals_utils.mod \
                gcxctbl_utils.mod kinds.mod lsd_func_utils.mod \
                metafun_utils.mod mp_interface.mod nvtx_utils.mod \
                parac.mod strs.mod system.mod tauf.mod tbxc.mod \
                timer.mod zeroing_utils.mod

gcxctbl_utils.mod.f90:$(SRCDIR)/gcxctbl_utils.mod.F90
gcxctbl_utils.mod.o:gcxctbl_utils.mod.f90 cnst.mod error_handling.mod \
                kinds.mod system.mod

genxc_utils.mod.f90:$(SRCDIR)/genxc_utils.mod.F90
genxc_utils.mod.o:genxc_utils.mod.f90 error_handling.mod functionals_utils.mod \
                kinds.mod parac.mod prmem_utils.mod tbxc.mod

geofile_utils.mod.f90:$(SRCDIR)/geofile_utils.mod.F90
geofile_utils.mod.o:geofile_utils.mod.f90 adat.mod cnst.mod \
                error_handling.mod fileopen_utils.mod fileopenmod.mod \
                ions.mod kinds.mod linres.mod meta_multiple_walkers_utils.mod \
                metr.mod mm_dim_utils.mod mm_dimmod.mod mw.mod \
                parac.mod pimd.mod readsr_utils.mod system.mod \
                bicanonicalCpmd.mod velupi_utils.mod

geq0.mod.f90:   $(SRCDIR)/geq0.mod.F90
geq0.mod.o:     geq0.mod.f90

getcor_utils.mod.f90:$(SRCDIR)/getcor_utils.mod.F90
getcor_utils.mod.o:getcor_utils.mod.f90 ions.mod kinds.mod \
                pbc_utils.mod pimd.mod system.mod zeroing_utils.mod

getfnm_utils.mod.f90:$(SRCDIR)/getfnm_utils.mod.F90
getfnm_utils.mod.o:getfnm_utils.mod.f90 cnst.mod ions.mod kinds.mod \
                pimd.mod system.mod zeroing_utils.mod

getfu_utils.mod.f90:$(SRCDIR)/getfu_utils.mod.F90
getfu_utils.mod.o:getfu_utils.mod.f90 cnst.mod ions.mod kinds.mod \
                pimd.mod system.mod

getgyr_utils.mod.f90:$(SRCDIR)/getgyr_utils.mod.F90
getgyr_utils.mod.o:getgyr_utils.mod.f90 cnst.mod ekinpp_utils.mod \
                ions.mod kinds.mod nose.mod pbc_utils.mod pimd.mod \
                system.mod zeroing_utils.mod

gettrans_utils.mod.f90:$(SRCDIR)/gettrans_utils.mod.F90
gettrans_utils.mod.o:gettrans_utils.mod.f90 cnst.mod dotp_utils.mod \
                error_handling.mod kinds.mod linres.mod mp_interface.mod \
                orbrot_utils.mod parac.mod spin.mod system.mod \
                zeroing_utils.mod

gfft_utils.mod.f90:$(SRCDIR)/gfft_utils.mod.F90
gfft_utils.mod.o:gfft_utils.mod.f90 error_handling.mod kinds.mod \
                parac.mod timer.mod kinds.mod error_handling.mod \
                timer.mod kinds.mod error_handling.mod timer.mod \
                kinds.mod error_handling.mod timer.mod

ghermit_utils.mod.f90:$(SRCDIR)/ghermit_utils.mod.F90
ghermit_utils.mod.o:ghermit_utils.mod.f90 error_handling.mod \
                kinds.mod

gle.mod.f90:    $(SRCDIR)/gle.mod.F90
gle.mod.o:      gle.mod.f90 kinds.mod

gle_utils.mod.f90:$(SRCDIR)/gle_utils.mod.F90
gle_utils.mod.o:gle_utils.mod.f90 cnst.mod coor.mod cotr.mod \
                error_handling.mod glemod.mod ions.mod isos.mod \
                kinds.mod mp_interface.mod parac.mod prng_utils.mod \
                rmas.mod store_types.mod system.mod zeroing_utils.mod

global_utils.mod.f90:$(SRCDIR)/global_utils.mod.F90
global_utils.mod.o:global_utils.mod.f90 kinds.mod mp_interface.mod \
                parac.mod pimd.mod

g_loc_dr_utils.mod.f90:$(SRCDIR)/g_loc_dr_utils.mod.F90
g_loc_dr_utils.mod.o:g_loc_dr_utils.mod.f90 cnst.mod cppt.mod \
                ddipo_utils.mod error_handling.mod fileopen_utils.mod \
                fileopenmod.mod g_loc.mod g_loc_exp_sum_utils.mod \
                g_loc_optim_utils.mod g_loc_realspace_utils.mod \
                g_loc_spread_ide_utils.mod g_loc_spread_sum_utils.mod \
                g_loc_util_utils.mod geq0mod.mod kinds.mod \
                kpnt.mod kpts.mod mp_interface.mod parac.mod \
                readsr_utils.mod setirec_utils.mod system.mod \
                utils.mod wann.mod wannier_center_utils.mod \
                wv30_utils.mod zeroing_utils.mod

g_loc_exp_sum_utils.mod.f90:$(SRCDIR)/g_loc_exp_sum_utils.mod.F90
g_loc_exp_sum_utils.mod.o:g_loc_exp_sum_utils.mod.f90 cnst.mod \
                cppt.mod ddipo_utils.mod error_handling.mod \
                fileopen_utils.mod fileopenmod.mod g_loc.mod \
                g_loc_optim_utils.mod g_loc_util_utils.mod \
                geq0mod.mod kinds.mod localize_utils.mod mp_interface.mod \
                parac.mod prng_utils.mod soft.mod system.mod \
                testex_utils.mod timer.mod u_upd_exp_sum_utils.mod \
                utils.mod wann.mod wannier_center_utils.mod \
                zeroing_utils.mod znum_mat_utils.mod

g_loc.mod.f90:  $(SRCDIR)/g_loc.mod.F90
g_loc.mod.o:    g_loc.mod.f90 kinds.mod

g_loc_opeigr_utils.mod.f90:$(SRCDIR)/g_loc_opeigr_utils.mod.F90
g_loc_opeigr_utils.mod.o:g_loc_opeigr_utils.mod.f90 kinds.mod \
                error_handling.mod timer.mod mp_interface.mod \
                opeigr_utils.mod system.mod parac.mod ddip.mod \
                g_loc.mod reshaper.mod zeroing_utils.mod

g_loc_optim_utils.mod.f90:$(SRCDIR)/g_loc_optim_utils.mod.F90
g_loc_optim_utils.mod.o:g_loc_optim_utils.mod.f90 cnst.mod \
                cppt.mod ddip.mod ddipo_utils.mod elct.mod \
                error_handling.mod fileopen_utils.mod fileopenmod.mod \
                g_loc.mod g_loc_util_utils.mod geq0mod.mod \
                kinds.mod localize_utils.mod mp_interface.mod \
                opeigr_utils.mod parac.mod prng_utils.mod soft.mod \
                system.mod testex_utils.mod timer.mod u_upd_exp_utils.mod \
                utils.mod wann.mod wannier_center_utils.mod \
                zeroing_utils.mod znum_mat_utils.mod kinds.mod \
                error_handling.mod timer.mod parac.mod kinds.mod \
                error_handling.mod timer.mod parac.mod g_loc.mod \
                kinds.mod error_handling.mod timer.mod system.mod \
                parac.mod g_loc.mod wann.mod

g_loc_realspace_utils.mod.f90:$(SRCDIR)/g_loc_realspace_utils.mod.F90
g_loc_realspace_utils.mod.o:g_loc_realspace_utils.mod.f90 coor.mod \
                cppt.mod error_handling.mod fftmain_utils.mod \
                fftutil_utils.mod fileopen_utils.mod fileopenmod.mod \
                g_loc.mod ions.mod kinds.mod mp_interface.mod \
                nmr_position_p_utils.mod parac.mod readsr_utils.mod \
                system.mod timer.mod zeroing_utils.mod

g_loc_spread_ide_utils.mod.f90:$(SRCDIR)/g_loc_spread_ide_utils.mod.F90
g_loc_spread_ide_utils.mod.o:g_loc_spread_ide_utils.mod.f90 \
                cnst.mod cppt.mod error_handling.mod fileopen_utils.mod \
                fileopenmod.mod g_loc.mod g_loc_util_utils.mod \
                geq0mod.mod kinds.mod metr.mod mp_interface.mod \
                parac.mod prng_utils.mod soft.mod system.mod \
                testex_utils.mod timer.mod utils.mod zeroing_utils.mod \
                kinds.mod error_handling.mod timer.mod system.mod \
                parac.mod kinds.mod error_handling.mod timer.mod \
                system.mod parac.mod zeroing_utils.mod

g_loc_spread_sum_utils.mod.f90:$(SRCDIR)/g_loc_spread_sum_utils.mod.F90
g_loc_spread_sum_utils.mod.o:g_loc_spread_sum_utils.mod.f90 \
                error_handling.mod fileopen_utils.mod fileopenmod.mod \
                g_loc.mod g_loc_spread_ide_utils.mod g_loc_util_utils.mod \
                geq0mod.mod kinds.mod mp_interface.mod parac.mod \
                prng_utils.mod soft.mod system.mod testex_utils.mod \
                timer.mod u_upd_spread_sum_utils.mod utils.mod \
                zeroing_utils.mod

g_loc_util_utils.mod.f90:$(SRCDIR)/g_loc_util_utils.mod.F90
g_loc_util_utils.mod.o:g_loc_util_utils.mod.f90 cnst.mod cppt.mod \
                error_handling.mod fileopen_utils.mod fileopenmod.mod \
                g_loc.mod geq0mod.mod kinds.mod mp_interface.mod \
                parac.mod prng_utils.mod rgs_utils.mod system.mod \
                timer.mod utils.mod wann.mod zeroing_utils.mod

glopar_utils.mod.f90:$(SRCDIR)/glopar_utils.mod.F90
glopar_utils.mod.o:glopar_utils.mod.f90 loadpa_utils.mod system.mod \
                timer.mod

gmopts_utils.mod.f90:$(SRCDIR)/gmopts_utils.mod.F90
gmopts_utils.mod.o:gmopts_utils.mod.f90 atwf.mod bsym.mod coor.mod \
                ddip.mod debfor_utils.mod elct.mod error_handling.mod \
                fnlalloc_utils.mod kinds.mod linres.mod lr_in_utils.mod \
                parac.mod prmem_utils.mod pslo.mod system.mod \
                timer.mod utils.mod vdwcmod.mod zeroing_utils.mod

gndstate_p_utils.mod.f90:$(SRCDIR)/gndstate_p_utils.mod.F90
gndstate_p_utils.mod.o:gndstate_p_utils.mod.f90 cppt.mod dotp_utils.mod \
                error_handling.mod fileopen_utils.mod fileopenmod.mod \
                gvec.mod kinds.mod mp_interface.mod nmr_position_p_utils.mod \
                parac.mod response_pmod.mod system.mod timer.mod \
                zeroing_utils.mod kinds.mod error_handling.mod \
                timer.mod system.mod parac.mod fileopenmod.mod \
                gndstate_p_utils.mod fileopen_utils.mod fileopen_utils.mod \
                fft_maxfft.mod

GPUBandwidth.f90:$(SRCDIR)/GPUBandwidth.F90
GPUBandwidth.o: GPUBandwidth.f90 cuda_types.mod cuda_utils.mod \
                kinds.mod sizeof_kinds.mod

graden_utils.mod.f90:$(SRCDIR)/graden_utils.mod.F90
graden_utils.mod.o:graden_utils.mod.f90 cnst.mod cppt.mod fft_maxfft.mod \
                fftmain_utils.mod fftnew_utils.mod kinds.mod \
                mp_interface.mod nvtx_utils.mod parac.mod state_utils.mod \
                system.mod timer.mod utils.mod zeroing_utils.mod

gs_disortho_utils.mod.f90:$(SRCDIR)/gs_disortho_utils.mod.F90
gs_disortho_utils.mod.o:gs_disortho_utils.mod.f90 error_handling.mod \
                geq0mod.mod kinds.mod mp_interface.mod parac.mod \
                system.mod timer.mod

gsize_utils.mod.f90:$(SRCDIR)/gsize_utils.mod.F90
gsize_utils.mod.o:gsize_utils.mod.f90 ions.mod kinds.mod

gsortho_utils.mod.f90:$(SRCDIR)/gsortho_utils.mod.F90
gsortho_utils.mod.o:gsortho_utils.mod.f90 dotp_utils.mod error_handling.mod \
                geq0mod.mod kinds.mod mp_interface.mod parac.mod \
                rgs_utils.mod system.mod timer.mod zeroing_utils.mod

gvec.mod.f90:   $(SRCDIR)/gvec.mod.F90
gvec.mod.o:     gvec.mod.f90 kinds.mod

h0psi1_p_utils.mod.f90:$(SRCDIR)/h0psi1_p_utils.mod.F90
h0psi1_p_utils.mod.o:h0psi1_p_utils.mod.f90 fft_maxfft.mod \
                fnonloc_utils.mod kinds.mod perturbation_p_utils.mod \
                response_pmod.mod rnlsm_utils.mod spin.mod \
                system.mod timer.mod vpsi_p_utils.mod

hardness_p_utils.mod.f90:$(SRCDIR)/hardness_p_utils.mod.F90
hardness_p_utils.mod.o:hardness_p_utils.mod.f90 coor.mod dotp_utils.mod \
                error_handling.mod fft_maxfft.mod fileopen_utils.mod \
                fileopenmod.mod geq0mod.mod ions.mod kinds.mod \
                mp_interface.mod parac.mod readsr_utils.mod \
                response_pmod.mod rhoofr_p_utils.mod rhoofr_utils.mod \
                rwfopt_p_utils.mod soft.mod spin.mod system.mod \
                zeroing_utils.mod

harm.mod.f90:   $(SRCDIR)/harm.mod.F90
harm.mod.o:     harm.mod.f90 kinds.mod

header_utils.mod.f90:$(SRCDIR)/header_utils.mod.F90
header_utils.mod.o:header_utils.mod.f90 envj.mod parac.mod \
                readsr_utils.mod

head.mod.f90:   $(SRCDIR)/head.mod.F90
head.mod.o:     head.mod.f90

hesele_p_utils.mod.f90:$(SRCDIR)/hesele_p_utils.mod.F90
hesele_p_utils.mod.o:hesele_p_utils.mod.f90 cppt.mod cvan.mod \
                ions.mod kinds.mod nlps.mod pslo.mod response_pmod.mod \
                sgpp.mod system.mod zeroing_utils.mod

hesele_utils.mod.f90:$(SRCDIR)/hesele_utils.mod.F90
hesele_utils.mod.o:hesele_utils.mod.f90 cppt.mod cvan.mod ions.mod \
                kinds.mod mp_interface.mod nlps.mod parac.mod \
                pslo.mod sgpp.mod simulmod.mod system.mod

hess_eta_p_utils.mod.f90:$(SRCDIR)/hess_eta_p_utils.mod.F90
hess_eta_p_utils.mod.o:hess_eta_p_utils.mod.f90 adat.mod d_mat_p_utils.mod \
                eicalc_utils.mod elct.mod error_handling.mod \
                fft_maxfft.mod fftmain_utils.mod fnonloc_p_utils.mod \
                implhv.mod ions.mod kinds.mod kpnt.mod nlps.mod \
                parac.mod response_pmod.mod rhoofr_p_utils.mod \
                rnlsm_p_utils.mod rnlsm_utils.mod rwfopt_p_utils.mod \
                sfac.mod spin.mod symm.mod system.mod timer.mod \
                zeroing_utils.mod

hessin_utils.mod.f90:$(SRCDIR)/hessin_utils.mod.F90
hessin_utils.mod.o:hessin_utils.mod.f90 coor.mod cotr.mod empfor_utils.mod \
                error_handling.mod fileopen_utils.mod fileopenmod.mod \
                ions.mod kinds.mod lscal.mod mp_interface.mod \
                parac.mod readsr_utils.mod symtrz_utils.mod \
                system.mod timer.mod utils.mod

hessout_utils.mod.f90:$(SRCDIR)/hessout_utils.mod.F90
hessout_utils.mod.o:hessout_utils.mod.f90 cotr.mod fileopen_utils.mod \
                fileopenmod.mod parac.mod

hessup_utils.mod.f90:$(SRCDIR)/hessup_utils.mod.F90
hessup_utils.mod.o:hessup_utils.mod.f90 cotr.mod error_handling.mod \
                kinds.mod timer.mod

hfx_drivers.mod.f90:$(SRCDIR)/hfx_drivers.mod.F90
hfx_drivers.mod.o:hfx_drivers.mod.f90 hfx_utils.mod hfxmod.mod \
                kinds.mod pw_hfx.mod pw_hfx_resp.mod pw_hfx_resp_types.mod \
                timer.mod

hfx.mod.f90:    $(SRCDIR)/hfx.mod.F90
hfx.mod.o:      hfx.mod.f90 kinds.mod

hfx_utils.mod.f90:$(SRCDIR)/hfx_utils.mod.F90
hfx_utils.mod.o:hfx_utils.mod.f90 cnst.mod cppt.mod dotp_utils.mod \
                error_handling.mod fft.mod fft_maxfft.mod fftmain_utils.mod \
                fftnew_utils.mod func.mod geq0mod.mod hfxmod.mod \
                kinds.mod kpts.mod machine.mod mp_interface.mod \
                parac.mod part_1d.mod pbc_utils.mod pslo.mod \
                ropt.mod rswfmod.mod spin.mod state_utils.mod \
                system.mod timer.mod zeroing_utils.mod

hipin_utils.mod.f90:$(SRCDIR)/hipin_utils.mod.F90
hipin_utils.mod.o:hipin_utils.mod.f90 cnst.mod cppt.mod error_handling.mod \
                fftmain_utils.mod geq0mod.mod gvec.mod isos.mod \
                kinds.mod parac.mod prmem_utils.mod special_functions.mod \
                system.mod timer.mod zeroing_utils.mod

hip_utils.mod.f90:$(SRCDIR)/hip_utils.mod.F90
hip_utils.mod.o:hip_utils.mod.f90 cppt.mod error_handling.mod \
                fft.mod fftmain_utils.mod isos.mod kinds.mod \
                parac.mod system.mod timer.mod zeroing_utils.mod

hnlmat_utils.mod.f90:$(SRCDIR)/hnlmat_utils.mod.F90
hnlmat_utils.mod.o:hnlmat_utils.mod.f90 cvan.mod error_handling.mod \
                ions.mod kinds.mod nlps.mod parac.mod pslo.mod \
                sfac.mod sgpp.mod spin.mod system.mod timer.mod

hpsi_utils.mod.f90:$(SRCDIR)/hpsi_utils.mod.F90
hpsi_utils.mod.o:hpsi_utils.mod.f90 coor.mod error_handling.mod \
                hubbardu.mod hubbardu_utils.mod fft_maxfft.mod \
                fnonloc_utils.mod geq0mod.mod kinds.mod kpts.mod \
                pslo.mod rnlsm_utils.mod spin.mod spsi_utils.mod \
                system.mod timer.mod utils.mod vpsi_utils.mod \
                zeroing_utils.mod

htrstr_utils.mod.f90:$(SRCDIR)/htrstr_utils.mod.F90
htrstr_utils.mod.o:htrstr_utils.mod.f90 cnst.mod cppt.mod geq0mod.mod \
                ions.mod kinds.mod pslo.mod ragg.mod sfac.mod \
                str2.mod strs.mod system.mod timer.mod

hubbardu.mod.f90:$(SRCDIR)/hubbardu.mod.F90
hubbardu.mod.o: hubbardu.mod.f90 kinds.mod

hubbardu_utils.mod.f90:$(SRCDIR)/hubbardu_utils.mod.F90
hubbardu_utils.mod.o:hubbardu_utils.mod.f90 atom.mod atwf.mod \
                cnst.mod cppt.mod cvan.mod hubbardu.mod domdr_utils.mod \
                elct.mod error_handling.mod fileopenmod.mod \
                fileopen_utils.mod fitpack_utils.mod fnlalloc_utils.mod \
                geq0mod.mod gvec.mod inscan_utils.mod ions.mod \
                kinds.mod kpts.mod lsfbtr_utils.mod mm_dimmod.mod \
                mp_interface.mod nlps.mod ortho_utils.mod ovlap_utils.mod \
                parac.mod phfac_utils.mod pslo.mod qspl.mod \
                readsr_utils.mod recpnew_utils.mod rnlsm_utils.mod \
                ropt.mod sfac.mod spin.mod spsi_utils.mod sphe.mod \
                system.mod timer.mod vdbp.mod ylmr_utils.mod \
                zeroing_utils.mod

if_parallel.mod.f90:$(SRCDIR)/if_parallel.mod.F90
if_parallel.mod.o:if_parallel.mod.f90

implhv.mod.f90: $(SRCDIR)/implhv.mod.F90
implhv.mod.o:   implhv.mod.f90 kinds.mod

initclust_utils.mod.f90:$(SRCDIR)/initclust_utils.mod.F90
initclust_utils.mod.o:initclust_utils.mod.f90 cnst.mod cppt.mod \
                cp_xc_utils.mod elct.mod error_handling.mod \
                fft.mod fftnew_utils.mod func.mod geq0mod.mod \
                hfxmod.mod hipin_utils.mod isos.mod kinds.mod \
                kpts.mod mp_interface.mod mtin_utils.mod parac.mod \
                rswfmod.mod spin.mod system.mod timer.mod

initrun_driver.mod.f90:$(SRCDIR)/initrun_driver.mod.F90
initrun_driver.mod.o:initrun_driver.mod.f90 ainitwf_utils.mod \
                andp.mod andr.mod bsym.mod calc_alm_utils.mod \
                chksym_utils.mod coor.mod copot_utils.mod elct.mod \
                error_handling.mod fint.mod geofile_utils.mod \
                isos.mod kinds.mod kpts.mod linres.mod merge_utils.mod \
                metr.mod mm_dim_utils.mod mm_dimmod.mod mm_input.mod \
                mp_interface.mod newcell_utils.mod nlcc.mod \
                ortho_utils.mod parac.mod phfac_utils.mod pimd.mod \
                pslo.mod ranc_utils.mod ranp_utils.mod rhoofr_c_utils.mod \
                rhoofr_utils.mod rinitwf_driver.mod rnlsm_utils.mod \
                ropt.mod rrandd_utils.mod rrane_utils.mod rv30_utils.mod \
                setirec_utils.mod setsc_utils.mod shop.mod \
                spin.mod store_types.mod system.mod timer.mod \
                vdw_utils.mod vdwcmod.mod wrgeo_utils.mod zeroing_utils.mod

initrun_utils.mod.f90:$(SRCDIR)/initrun_utils.mod.F90
initrun_utils.mod.o:initrun_utils.mod.f90 ainitwf_utils.mod \
                calc_alm_utils.mod copot_utils.mod elct.mod \
                fint.mod newcell_utils.mod nlcc.mod ortho_utils.mod \
                pslo.mod rhoofr_utils.mod rinitwf_utils.mod \
                rnlsm_utils.mod store_types.mod system.mod \
                timer.mod

inr_dr_utils.mod.f90:$(SRCDIR)/inr_dr_utils.mod.F90
inr_dr_utils.mod.o:inr_dr_utils.mod.f90 coor.mod cotr.mod do_perturbation_p_utils.mod \
                eicalc_utils.mod error_handling.mod fft_maxfft.mod \
                fileopen_utils.mod fileopenmod.mod geofile_utils.mod \
                gsize_utils.mod hess_eta_p_utils.mod hessup_utils.mod \
                implhv.mod ions.mod kinds.mod lanc_phon_p_utils.mod \
                mp_interface.mod nlps.mod norm.mod parac.mod \
                perturbation_p_utils.mod phonons_p_utils.mod \
                puttau_utils.mod response_pmod.mod rgdiis_utils.mod \
                rhoofr_utils.mod rmas.mod rnlsm_p_utils.mod \
                rnlsm_utils.mod ropt.mod sdion_utils.mod sfac.mod \
                soft.mod system.mod testex_utils.mod timer.mod \
                utils.mod vofrho_utils.mod wrgeo_utils.mod \
                xinr.mod zeroing_utils.mod

inscan_utils.mod.f90:$(SRCDIR)/inscan_utils.mod.F90
inscan_utils.mod.o:inscan_utils.mod.f90 error_handling.mod \
                readsr_utils.mod store_types.mod

interaction_manno_p_utils.mod.f90:$(SRCDIR)/interaction_manno_p_utils.mod.F90
interaction_manno_p_utils.mod.o:interaction_manno_p_utils.mod.f90 \
                atomwf_utils.mod atwf.mod coor.mod cppt.mod \
                csize_utils.mod ddip.mod ddipo_utils.mod eicalc_utils.mod \
                elct.mod ener.mod error_handling.mod fft_maxfft.mod \
                fileopen_utils.mod fileopenmod.mod forces_driver.mod \
                forces_utils.mod kinds.mod localize_utils.mod \
                lowdin_utils.mod machine.mod opeigr_utils.mod \
                ovlap_utils.mod parac.mod perturbation_p_utils.mod \
                prmem_utils.mod prop.mod readsr_utils.mod response_pmod.mod \
                rhoofr_utils.mod ropt.mod rotate_my_wannier_manno_p_utils.mod \
                rwfopt_p_utils.mod setirec_utils.mod store_types.mod \
                system.mod timer.mod utils.mod vofrho_utils.mod \
                wann.mod wannier_center_utils.mod wv30_utils.mod \
                zeroing_utils.mod kinds.mod error_handling.mod \
                timer.mod system.mod parac.mod coor.mod fft_maxfft.mod \
                zeroing_utils.mod

interaction_p_utils.mod.f90:$(SRCDIR)/interaction_p_utils.mod.F90
interaction_p_utils.mod.o:interaction_p_utils.mod.f90 atomwf_utils.mod \
                atwf.mod coor.mod csize_utils.mod ddip.mod \
                eicalc_utils.mod ener.mod error_handling.mod \
                fft_maxfft.mod forces_driver.mod kinds.mod \
                localize_utils.mod lowdin_utils.mod machine.mod \
                parac.mod perturbation_p_utils.mod prmem_utils.mod \
                prop.mod readsr_utils.mod response_pmod.mod \
                rhoofr_utils.mod ropt.mod rotate_my_wannier_para_p_utils.mod \
                rwfopt_p_utils.mod sdlinres_utils.mod setirec_utils.mod \
                sort_utils.mod store_types.mod system.mod timer.mod \
                utils.mod vofrho_utils.mod wann.mod wv30_utils.mod \
                zeroing_utils.mod kinds.mod error_handling.mod \
                timer.mod mp_interface.mod system.mod parac.mod \
                coor.mod geq0mod.mod response_pmod.mod fft_maxfft.mod \
                readsr_utils.mod interaction_p_utils.mod zeroing_utils.mod

interp3d_utils.mod.f90:$(SRCDIR)/interp3d_utils.mod.F90
interp3d_utils.mod.o:interp3d_utils.mod.f90 error_handling.mod \
                kinds.mod parac.mod zeroing_utils.mod

interpt_utils.mod.f90:$(SRCDIR)/interpt_utils.mod.F90
interpt_utils.mod.o:interpt_utils.mod.f90 ddip.mod egointer_utils.mod \
                elct.mod error_handling.mod fnlalloc_utils.mod \
                kinds.mod parac.mod prmem_utils.mod pslo.mod \
                rwfopt_utils.mod system.mod utils.mod

ions.mod.f90:   $(SRCDIR)/ions.mod.F90
ions.mod.o:     ions.mod.f90 kinds.mod system.mod

io_utils.mod.f90:$(SRCDIR)/io_utils.mod.F90
io_utils.mod.o: io_utils.mod.f90 error_handling.mod kinds.mod \
                string_utils.mod

isos.mod.f90:   $(SRCDIR)/isos.mod.F90
isos.mod.o:     isos.mod.f90 kinds.mod

jacobi_c_utils.mod.f90:$(SRCDIR)/jacobi_c_utils.mod.F90
jacobi_c_utils.mod.o:jacobi_c_utils.mod.f90 kinds.mod timer.mod

jacobi_utils.mod.f90:$(SRCDIR)/jacobi_utils.mod.F90
jacobi_utils.mod.o:jacobi_utils.mod.f90 kinds.mod timer.mod

jrotation_utils.mod.f90:$(SRCDIR)/jrotation_utils.mod.F90
jrotation_utils.mod.o:jrotation_utils.mod.f90 cnst.mod error_handling.mod \
                kinds.mod machine.mod mp_interface.mod parac.mod \
                prng_utils.mod sd_wannier_utils.mod system.mod \
                timer.mod utils.mod wann.mod zeroing_utils.mod

k290_2_utils.mod.f90:$(SRCDIR)/k290_2_utils.mod.F90
k290_2_utils.mod.o:k290_2_utils.mod.f90 error_handling.mod \
                kinds.mod parac.mod readsr_utils.mod utils.mod

k290_utils.mod.f90:$(SRCDIR)/k290_utils.mod.F90
k290_utils.mod.o:k290_utils.mod.f90 error_handling.mod k290_2_utils.mod \
                kinds.mod parac.mod

kddipo_utils.mod.f90:$(SRCDIR)/kddipo_utils.mod.F90
kddipo_utils.mod.o:kddipo_utils.mod.f90 ddipo_utils.mod error_handling.mod \
                kinds.mod kpnt.mod kpts.mod opeigr_c_utils.mod \
                opeigr_utils.mod parac.mod rggen_utils.mod \
                system.mod timer.mod

k_diis_rhofix_utils.mod.f90:$(SRCDIR)/k_diis_rhofix_utils.mod.F90
k_diis_rhofix_utils.mod.o:k_diis_rhofix_utils.mod.f90 elct.mod \
                ener.mod error_handling.mod k_forces_driver.mod \
                kinds.mod kpnt.mod machine.mod norm.mod parac.mod \
                ropt.mod soft.mod store_types.mod system.mod \
                testex_utils.mod wrener_utils.mod zeroing_utils.mod

kdpc.mod.f90:   $(SRCDIR)/kdpc.mod.F90
kdpc.mod.o:     kdpc.mod.f90 kinds.mod

kdp_diag_utils.mod.f90:$(SRCDIR)/kdp_diag_utils.mod.F90
kdp_diag_utils.mod.o:kdp_diag_utils.mod.f90 cnst.mod error_handling.mod \
                kinds.mod parac.mod ropt.mod system.mod timer.mod

kdp.mod.f90:    $(SRCDIR)/kdp.mod.F90
kdp.mod.o:      kdp.mod.f90 kinds.mod

kdpoints_utils.mod.f90:$(SRCDIR)/kdpoints_utils.mod.F90
kdpoints_utils.mod.o:kdpoints_utils.mod.f90 error_handling.mod \
                kdp.mod kdpc.mod kpnt.mod mp_interface.mod \
                parac.mod rkpnt_utils.mod system.mod

kdp_prep_utils.mod.f90:$(SRCDIR)/kdp_prep_utils.mod.F90
kdp_prep_utils.mod.o:kdp_prep_utils.mod.f90 calc_pij_utils.mod \
                error_handling.mod kinds.mod mp_interface.mod \
                parac.mod spin.mod system.mod zeroing_utils.mod

kdp_rho_utils.mod.f90:$(SRCDIR)/kdp_rho_utils.mod.F90
kdp_rho_utils.mod.o:kdp_rho_utils.mod.f90 error_handling.mod \
                fft_maxfft.mod kinds.mod rhoofr_kdp_utils.mod \
                spin.mod system.mod timer.mod zeroing_utils.mod

kdp_stress_kin_utils.mod.f90:$(SRCDIR)/kdp_stress_kin_utils.mod.F90
kdp_stress_kin_utils.mod.o:kdp_stress_kin_utils.mod.f90 error_handling.mod \
                kinds.mod parac.mod strs.mod system.mod zeroing_utils.mod

k_forces_driver.mod.f90:$(SRCDIR)/k_forces_driver.mod.F90
k_forces_driver.mod.o:k_forces_driver.mod.f90 andp.mod cnst.mod \
                csize_utils.mod dotp_utils.mod elct.mod ener.mod \
                error_handling.mod fnonloc_utils.mod forces_driver.mod \
                frsblk_c_utils.mod frsblk_utils.mod geq0mod.mod \
                gsize_utils.mod hesele_utils.mod hnlmat_utils.mod \
                k_hesele_utils.mod k_pcgrad_utils.mod kinds.mod \
                kpnt.mod kpts.mod mp_interface.mod nlforce_utils.mod \
                nlps.mod norm.mod ortho_utils.mod ovlap_utils.mod \
                parac.mod pcgrad_driver.mod pslo.mod puttau_utils.mod \
                reshaper.mod rhoofr_c_utils.mod rnlfl_utils.mod \
                rnlfor_utils.mod rnlrh_utils.mod rnlsm_utils.mod \
                ropt.mod rotate_utils.mod rpiiint_utils.mod \
                rswfmod.mod sfac.mod sort_utils.mod spin.mod \
                stress_utils.mod summat_utils.mod symtrz_utils.mod \
                system.mod timer.mod tpar.mod utils.mod vdw_utils.mod \
                vdwcmod.mod vofrho_utils.mod vpsi_utils.mod \
                zeroing_utils.mod

k_forces_utils.mod.f90:$(SRCDIR)/k_forces_utils.mod.F90
k_forces_utils.mod.o:k_forces_utils.mod.f90 fft_maxfft.mod \
                fnonloc_utils.mod nlforce_utils.mod nlps.mod \
                ortho_utils.mod pslo.mod rnlsm_utils.mod ropt.mod \
                rscpot_utils.mod spin.mod summat_utils.mod \
                symtrz_utils.mod system.mod

k_hesele_utils.mod.f90:$(SRCDIR)/k_hesele_utils.mod.F90
k_hesele_utils.mod.o:k_hesele_utils.mod.f90 cppt.mod ions.mod \
                kinds.mod kpnt.mod mp_interface.mod nlps.mod \
                parac.mod sgpp.mod simulmod.mod system.mod

kinds.mod.f90:  $(SRCDIR)/kinds.mod.F90
kinds.mod.o:    kinds.mod.f90

kin_energy_utils.mod.f90:$(SRCDIR)/kin_energy_utils.mod.F90
kin_energy_utils.mod.o:kin_energy_utils.mod.f90 cppt.mod dotp_utils.mod \
                elct.mod ener.mod error_handling.mod geq0mod.mod \
                kinds.mod prcp.mod special_functions.mod spin.mod \
                system.mod timer.mod

k_odiis_utils.mod.f90:$(SRCDIR)/k_odiis_utils.mod.F90
k_odiis_utils.mod.o:k_odiis_utils.mod.f90 zeroing_utils.mod \
                kinds.mod error_handling.mod timer.mod mp_interface.mod \
                system.mod parac.mod ener.mod elct.mod geq0mod.mod \
                odiis_utils.mod zeroing_utils.mod kinds.mod \
                error_handling.mod timer.mod mp_interface.mod \
                system.mod parac.mod

k_pcgrad_utils.mod.f90:$(SRCDIR)/k_pcgrad_utils.mod.F90
k_pcgrad_utils.mod.o:k_pcgrad_utils.mod.f90 dotp_utils.mod \
                elct.mod ener.mod error_handling.mod kinds.mod \
                mp_interface.mod ortho_utils.mod parac.mod \
                pcgrad_driver.mod pslo.mod rnlsm_utils.mod \
                system.mod timer.mod tpar.mod zeroing_utils.mod

kpclean_utils.mod.f90:$(SRCDIR)/kpclean_utils.mod.F90
kpclean_utils.mod.o:kpclean_utils.mod.f90 kinds.mod sphe.mod \
                system.mod

kpert_potential_p_utils.mod.f90:$(SRCDIR)/kpert_potential_p_utils.mod.F90
kpert_potential_p_utils.mod.o:kpert_potential_p_utils.mod.f90 \
                cppt.mod error_handling.mod ions.mod kinds.mod \
                nlps.mod parac.mod pslo.mod response_pmod.mod \
                sfac.mod sgpp.mod system.mod timer.mod zeroing_utils.mod

kpert_util_p_utils.mod.f90:$(SRCDIR)/kpert_util_p_utils.mod.F90
kpert_util_p_utils.mod.o:kpert_util_p_utils.mod.f90 csmat_utils.mod \
                error_handling.mod kinds.mod kpnt.mod response_pmod.mod \
                sfac.mod system.mod utils.mod zeroing_utils.mod

kpnt.mod.f90:   $(SRCDIR)/kpnt.mod.F90
kpnt.mod.o:     kpnt.mod.f90 kinds.mod

kpts.mod.f90:   $(SRCDIR)/kpts.mod.F90
kpts.mod.o:     kpts.mod.f90 kinds.mod

ksdiag_utils.mod.f90:$(SRCDIR)/ksdiag_utils.mod.F90
ksdiag_utils.mod.o:ksdiag_utils.mod.f90 cppt.mod cvan.mod ions.mod \
                kinds.mod mp_interface.mod nlps.mod parac.mod \
                pslo.mod sgpp.mod simulmod.mod system.mod

ks_ener_p_utils.mod.f90:$(SRCDIR)/ks_ener_p_utils.mod.F90
ks_ener_p_utils.mod.o:ks_ener_p_utils.mod.f90 cnst.mod error_handling.mod \
                frsblk_c_utils.mod kinds.mod kpert_util_p_utils.mod \
                kpnt.mod parac.mod response_pmod.mod system.mod \
                timer.mod zeroing_utils.mod

ksmat_dist_utils.mod.f90:$(SRCDIR)/ksmat_dist_utils.mod.F90
ksmat_dist_utils.mod.o:ksmat_dist_utils.mod.f90 atwf.mod error_handling.mod \
                fnlalloc_utils.mod fnonloc_utils.mod gsortho_utils.mod \
                ions.mod jrotation_utils.mod kinds.mod kpts.mod \
                ksmatmod.mod mp_interface.mod parac.mod prng_utils.mod \
                rnlsm_utils.mod spin.mod system.mod timer.mod \
                vpsi_utils.mod wfnio_utils.mod zeroing_utils.mod

ksmat.mod.f90:  $(SRCDIR)/ksmat.mod.F90
ksmat.mod.o:    ksmat.mod.f90 kinds.mod

ksmat_utils.mod.f90:$(SRCDIR)/ksmat_utils.mod.F90
ksmat_utils.mod.o:ksmat_utils.mod.f90 atwf.mod error_handling.mod \
                fnlalloc_utils.mod fnonloc_utils.mod ions.mod \
                kinds.mod kpts.mod rnlsm_utils.mod spin.mod \
                system.mod timer.mod vpsi_utils.mod zeroing_utils.mod

k_updwf_utils.mod.f90:$(SRCDIR)/k_updwf_utils.mod.F90
k_updwf_utils.mod.o:k_updwf_utils.mod.f90 elct.mod error_handling.mod \
                forcedr_utils.mod k_forces_driver.mod kinds.mod \
                kpnt.mod norm.mod ortho_utils.mod parac.mod \
                pcgrad_utils.mod ropt.mod soft.mod system.mod \
                testex_utils.mod zeroing_utils.mod

lanc_phon_p_utils.mod.f90:$(SRCDIR)/lanc_phon_p_utils.mod.F90
lanc_phon_p_utils.mod.o:lanc_phon_p_utils.mod.f90 coor.mod \
                cotr.mod d_mat_p_utils.mod eicalc_utils.mod \
                elct.mod error_handling.mod fft_maxfft.mod \
                fileopen_utils.mod fileopenmod.mod hess_eta_p_utils.mod \
                implhv.mod ions.mod kinds.mod kpnt.mod mp_interface.mod \
                nlps.mod parac.mod phonons_p_utils.mod prng_utils.mod \
                readsr_utils.mod response_pmod.mod rmas.mod \
                rnlsm_p_utils.mod rnlsm_utils.mod sfac.mod \
                soft.mod softex_utils.mod system.mod testex_utils.mod \
                timer.mod zeroing_utils.mod

latgen_utils.mod.f90:$(SRCDIR)/latgen_utils.mod.F90
latgen_utils.mod.o:latgen_utils.mod.f90 error_handling.mod \
                kinds.mod parac.mod

ldos.mod.f90:   $(SRCDIR)/ldos.mod.F90
ldos.mod.o:     ldos.mod.f90

ldos_utils.mod.f90:$(SRCDIR)/ldos_utils.mod.F90
ldos_utils.mod.o:ldos_utils.mod.f90 cnst.mod cppt.mod elct.mod \
                ener.mod error_handling.mod fft.mod fft_maxfft.mod \
                fftmain_utils.mod fileopen_utils.mod fileopenmod.mod \
                geq0mod.mod kinds.mod kpnt.mod kpts.mod ldosmod.mod \
                mp_interface.mod parac.mod system.mod timer.mod \
                zeroing_utils.mod

legendre_p_utils.mod.f90:$(SRCDIR)/legendre_p_utils.mod.F90
legendre_p_utils.mod.o:legendre_p_utils.mod.f90 error_handling.mod \
                kinds.mod parac.mod

linalg_utils.mod.f90:$(SRCDIR)/linalg_utils.mod.F90
linalg_utils.mod.o:linalg_utils.mod.f90 error_handling.mod \
                kinds.mod mp_interface.mod system.mod zeroing_utils.mod

linres.mod.f90: $(SRCDIR)/linres.mod.F90
linres.mod.o:   linres.mod.f90 kinds.mod

loadpa_utils.mod.f90:$(SRCDIR)/loadpa_utils.mod.F90
loadpa_utils.mod.o:loadpa_utils.mod.f90 cppt.mod elct.mod error_handling.mod \
                geq0mod.mod gvec.mod ions.mod isos.mod kinds.mod \
                kpts.mod mm_dim_utils.mod mm_dimmod.mod mp_interface.mod \
                parac.mod prmem_utils.mod sort_utils.mod sphe.mod \
                system.mod timer.mod zeroing_utils.mod

loadse_utils.mod.f90:$(SRCDIR)/loadse_utils.mod.F90
loadse_utils.mod.o:loadse_utils.mod.f90 cppt.mod elct.mod error_handling.mod \
                gvec.mod ions.mod kinds.mod mm_dim_utils.mod \
                mm_dimmod.mod parac.mod sphe.mod system.mod

localize_utils.mod.f90:$(SRCDIR)/localize_utils.mod.F90
localize_utils.mod.o:localize_utils.mod.f90 atwf.mod ddip.mod \
                ddipo_utils.mod dotp_utils.mod error_handling.mod \
                forcedr_driver.mod forcep_utils.mod gvec.mod \
                ions.mod jrotation_utils.mod kinds.mod linres.mod \
                molorb_utils.mod mp_interface.mod nmr_position_p_utils.mod \
                opeigr_utils.mod parac.mod response_pmod.mod \
                rmas.mod rotate_utils.mod setbasis_utils.mod \
                spin.mod system.mod timer.mod utils.mod wann.mod \
                wannier_center_utils.mod wannier_print_utils.mod \
                wc_dos_utils.mod zeroing_utils.mod

locpot.mod.f90: $(SRCDIR)/locpot.mod.F90
locpot.mod.o:   locpot.mod.f90 kinds.mod

lodipo_utils.mod.f90:$(SRCDIR)/lodipo_utils.mod.F90
lodipo_utils.mod.o:lodipo_utils.mod.f90 cnst.mod cppt.mod geq0mod.mod \
                kinds.mod lodp.mod mp_interface.mod parac.mod \
                system.mod timer.mod

lodp.mod.f90:   $(SRCDIR)/lodp.mod.F90
lodp.mod.o:     lodp.mod.f90 kinds.mod

lowdin_utils.mod.f90:$(SRCDIR)/lowdin_utils.mod.F90
lowdin_utils.mod.o:lowdin_utils.mod.f90 csmat_utils.mod error_handling.mod \
                geq0mod.mod kinds.mod parac.mod reshaper.mod \
                sfac.mod spin.mod system.mod timer.mod utils.mod \
                zeroing_utils.mod

lr_diag_utils.mod.f90:$(SRCDIR)/lr_diag_utils.mod.F90
lr_diag_utils.mod.o:lr_diag_utils.mod.f90 error_handling.mod \
                kinds.mod linres.mod system.mod td_dav_utils.mod \
                td_nhdav_utils.mod td_pcg_utils.mod timer.mod

lr_force_utils.mod.f90:$(SRCDIR)/lr_force_utils.mod.F90
lr_force_utils.mod.o:lr_force_utils.mod.f90 dotp_utils.mod \
                elct.mod error_handling.mod hfx_drivers.mod \
                hpsi_utils.mod kinds.mod lr_xcpot_utils.mod \
                mp_interface.mod ovlap_utils.mod parac.mod \
                poin.mod rho1ofr_utils.mod rotate_utils.mod \
                spin.mod system.mod timer.mod v1ofrho1_utils.mod \
                vpsi_utils.mod

lr_in_utils.mod.f90:$(SRCDIR)/lr_in_utils.mod.F90
lr_in_utils.mod.o:lr_in_utils.mod.f90 adat.mod coor.mod dftin_utils.mod \
                error_handling.mod inscan_utils.mod ions.mod \
                kinds.mod linres.mod lr_xcpot_utils.mod mols.mod \
                mp_interface.mod parac.mod pw_hfx_resp_types.mod \
                readsr_utils.mod soc_types.mod spin.mod system.mod \
                zeroing_utils.mod

lr_ortho_utils.mod.f90:$(SRCDIR)/lr_ortho_utils.mod.F90
lr_ortho_utils.mod.o:lr_ortho_utils.mod.f90 elct.mod error_handling.mod \
                geq0mod.mod kinds.mod linres.mod mp_interface.mod \
                ovlap_utils.mod parac.mod rotate_utils.mod \
                spin.mod system.mod timer.mod utils.mod

lr_pcg_utils.mod.f90:$(SRCDIR)/lr_pcg_utils.mod.F90
lr_pcg_utils.mod.o:lr_pcg_utils.mod.f90 cppt.mod dotp_utils.mod \
                elct.mod fftmain_utils.mod hpsi_utils.mod kinds.mod \
                linres.mod lr_ortho_utils.mod mp_interface.mod \
                parac.mod poin.mod rho1ofr_utils.mod spin.mod \
                system.mod timer.mod zeroing_utils.mod

lr_tddft_drhoe.mod.f90:$(SRCDIR)/lr_tddft_drhoe.mod.F90
lr_tddft_drhoe.mod.o:lr_tddft_drhoe.mod.f90 cppt.mod elct.mod \
                error_handling.mod fftmain_utils.mod kinds.mod \
                linres.mod mp_interface.mod ovlap_utils.mod \
                parac.mod pimd.mod poin.mod readsr_utils.mod \
                rho1ofr_utils.mod rotate_utils.mod spin.mod \
                system.mod

lr_tddft_utils.mod.f90:$(SRCDIR)/lr_tddft_utils.mod.F90
lr_tddft_utils.mod.o:lr_tddft_utils.mod.f90 canon_utils.mod \
                corec_utils.mod ddipo_utils.mod dftin_utils.mod \
                elct.mod ener.mod error_handling.mod forcedr_driver.mod \
                forcedr_utils.mod gettrans_utils.mod ions.mod \
                isos.mod kinds.mod linres.mod lr_diag_utils.mod \
                lr_tddft_drhoe.mod lr_xcpot_utils.mod mm_dim_utils.mod \
                mm_dimmod.mod mm_input.mod mm_qmmm_forcedr_utils.mod \
                mm_rho_forcedr_utils.mod mp_interface.mod nlcc.mod \
                ortho_utils.mod ovlap_utils.mod parac.mod poin.mod \
                rhoofr_utils.mod rnlsm_utils.mod ropt.mod sh_tddft_utils.mod \
                sh_utils.mod soc_types.mod spin.mod system.mod \
                tauf.mod td_force_utils.mod td_mm_qmmm_forcedr_utils.mod \
                td_os_utils.mod timer.mod tpot.mod zeroing_utils.mod

lr_upd_utils.mod.f90:$(SRCDIR)/lr_upd_utils.mod.F90
lr_upd_utils.mod.o:lr_upd_utils.mod.f90 csize_utils.mod dotp_utils.mod \
                elct.mod error_handling.mod geq0mod.mod kinds.mod \
                linres.mod lr_force_utils.mod lr_ortho_utils.mod \
                lr_pcg_utils.mod mp_interface.mod norm.mod \
                parac.mod ropt.mod soft.mod system.mod testex_utils.mod \
                timer.mod utils.mod zdiis_utils.mod

lr_xcpot_utils.mod.f90:$(SRCDIR)/lr_xcpot_utils.mod.F90
lr_xcpot_utils.mod.o:lr_xcpot_utils.mod.f90 cp_xc_utils.mod \
                density_functionals_utils.mod error_handling.mod \
                func.mod kinds.mod linres.mod nlcc.mod parac.mod \
                spin.mod system.mod tbxc.mod zeroing_utils.mod

lscal.mod.f90:  $(SRCDIR)/lscal.mod.F90
lscal.mod.o:    lscal.mod.f90 kinds.mod nvarmod.mod

lsd_elf_utils.mod.f90:$(SRCDIR)/lsd_elf_utils.mod.F90
lsd_elf_utils.mod.o:lsd_elf_utils.mod.f90 bsym.mod cnst.mod \
                cppt.mod elct.mod error_handling.mod fft_maxfft.mod \
                fftmain_utils.mod fftnew_utils.mod kinds.mod \
                kpts.mod meta_multiple_walkers_utils.mod mp_interface.mod \
                mw.mod parac.mod phfac_utils.mod pimd.mod prden.mod \
                pslo.mod readsr_utils.mod rhoofr_c_utils.mod \
                rhoofr_utils.mod rnlsm_utils.mod spin.mod system.mod \
                zeroing_utils.mod

lsd_func_utils.mod.f90:$(SRCDIR)/lsd_func_utils.mod.F90
lsd_func_utils.mod.o:lsd_func_utils.mod.f90 cnst.mod error_handling.mod \
                func.mod functionals_utils.mod kinds.mod system.mod

lsfbtr_utils.mod.f90:$(SRCDIR)/lsfbtr_utils.mod.F90
lsfbtr_utils.mod.o:lsfbtr_utils.mod.f90 kinds.mod

lsforce_utils.mod.f90:$(SRCDIR)/lsforce_utils.mod.F90
lsforce_utils.mod.o:lsforce_utils.mod.f90 bsym.mod kinds.mod \
                system.mod

lxc_utils.mod.f90:$(SRCDIR)/lxc_utils.mod.F90
lxc_utils.mod.o:lxc_utils.mod.f90 error_handling.mod kinds.mod

machine.mod.f90:$(SRCDIR)/machine.mod.F90
machine.mod.o:  machine.mod.f90 kinds.mod

matrix_p_utils.mod.f90:$(SRCDIR)/matrix_p_utils.mod.F90
matrix_p_utils.mod.o:matrix_p_utils.mod.f90 cppt.mod error_handling.mod \
                fft_maxfft.mod fnonloc_utils.mod geq0mod.mod \
                ions.mod kinds.mod machine.mod mp_interface.mod \
                nlps.mod ovlap_utils.mod parac.mod pslo.mod \
                response_pmod.mod rnlsm_utils.mod sfac.mod \
                sgpp.mod spin.mod system.mod timer.mod utils.mod \
                vpsi_utils.mod zeroing_utils.mod

mdclas_utils.mod.f90:$(SRCDIR)/mdclas_utils.mod.F90
mdclas_utils.mod.o:mdclas_utils.mod.f90 anneal_utils.mod clas.mod \
                clas_force_utils.mod cnst.mod comvel_utils.mod \
                comvelmod.mod coor.mod copot_utils.mod cotr.mod \
                csize_utils.mod ddipo_utils.mod deort_utils.mod \
                detdof_utils.mod dispp_utils.mod dynit_utils.mod \
                elct.mod ener.mod error_handling.mod fileopen_utils.mod \
                fileopenmod.mod finalp_utils.mod forcedr_driver.mod \
                forcep_utils.mod freqs_utils.mod geofile_utils.mod \
                geq0mod.mod gsize_utils.mod initrun_driver.mod \
                ions.mod kinds.mod machine.mod mdmain_utils.mod \
                metr.mod movi.mod mp_interface.mod nlcc.mod \
                norm.mod noseinit_utils.mod noseng_utils.mod \
                nosepa_utils.mod noseup_utils.mod nospinit_utils.mod \
                parac.mod phfac_utils.mod posupa_utils.mod \
                printave_utils.mod printp_utils.mod prng_utils.mod \
                proja_utils.mod pslo.mod quenbo_utils.mod rattle_utils.mod \
                rekine_utils.mod resetac_utils.mod rhopri_utils.mod \
                ropt.mod rortv_utils.mod rotvel_utils.mod rscve_utils.mod \
                setirec_utils.mod shake_utils.mod soft.mod \
                spin.mod store_types.mod system.mod temps.mod \
                testex_utils.mod totstr_utils.mod tpar.mod \
                utils.mod velupa_utils.mod wannier_print_utils.mod \
                wrener_utils.mod wrgeo_utils.mod wv30_utils.mod \
                zeroing_utils.mod

mddiag_interaction_p_utils.mod.f90:$(SRCDIR)/mddiag_interaction_p_utils.mod.F90
mddiag_interaction_p_utils.mod.o:mddiag_interaction_p_utils.mod.f90 \
                cppt.mod ddip.mod do_perturbation_p_utils.mod \
                eicalc_utils.mod ener.mod error_handling.mod \
                forcep_utils.mod forces_driver.mod kinds.mod \
                localize_utils.mod lowdin_utils.mod mp_interface.mod \
                parac.mod perturbation_p_utils.mod phfac_utils.mod \
                prmem_utils.mod prop.mod response_pmod.mod \
                rhoofr_utils.mod rwfopt_p_utils.mod spin.mod \
                system.mod utils.mod vofrho_utils.mod wann.mod \
                zeroing_utils.mod

md_driver.mod.f90:$(SRCDIR)/md_driver.mod.F90
md_driver.mod.o:md_driver.mod.f90 andp.mod andr.mod anneal_utils.mod \
                atwf.mod box_boundary_utils.mod bs_forces_diag_utils.mod \
                bsym.mod calc_alm_utils.mod cdft_utils.mod \
                cdftmod.mod cnst.mod cnst_dyn.mod comvel_utils.mod \
                comvelmod.mod coor.mod copot_utils.mod cotr.mod \
                ddipo_utils.mod detdof_utils.mod dispp_utils.mod \
                dynit_utils.mod efld.mod ekinpp_utils.mod elct.mod \
                ener.mod error_handling.mod extrap_utils.mod \
                fileopen_utils.mod fileopenmod.mod finalp_utils.mod \
                fint.mod forcep_utils.mod forces_diag_utils.mod \
                forces_prop_utils.mod geofile_utils.mod gle_utils.mod \
                glemod.mod gsize_utils.mod hfxmod.mod hubbardu.mod \
                initrun_driver.mod initrun_utils.mod ions.mod \
                isos.mod kinds.mod kpts.mod linres.mod localize_utils.mod \
                lr_tddft_utils.mod machine.mod mddiag_interaction_p_utils.mod \
                meta_colvar_inp_utils.mod meta_colvar_utils.mod \
                meta_exl_mult_utils.mod meta_exlagr_methods.mod \
                meta_exlagr_utils.mod meta_multiple_walkers_utils.mod \
                mfep.mod mm_extrap.mod moverho_utils.mod mp_interface.mod \
                mw.mod nlcc.mod norm.mod nose.mod noseng_utils.mod \
                nosepa_utils.mod noseup_utils.mod nospinit_utils.mod \
                ovlap_utils.mod parac.mod phfac_utils.mod pimd.mod \
                poin.mod posupi_utils.mod printave_utils.mod \
                printp_utils.mod proppt_utils.mod pslo.mod \
                puttau_utils.mod rattle_utils.mod readsr_utils.mod \
                resetac_utils.mod response_pmod.mod rhoofr_utils.mod \
                rhopri_utils.mod rinitwf_utils.mod rinvel_utils.mod \
                rnlsm_utils.mod ropt.mod rotate_utils.mod rotvel_utils.mod \
                rscvp_utils.mod sample_utils.mod setbsstate_utils.mod \
                setirec_utils.mod sh_tddft_utils.mod shake_utils.mod \
                soc.mod soc_types.mod soft.mod spin.mod store_types.mod \
                system.mod td_utils.mod testex_utils.mod teststore_utils.mod \
                timer.mod totstr_utils.mod vdwcmod.mod velupi_utils.mod \
                wrener_utils.mod wrgeo_utils.mod wv30_utils.mod \
                zeroing_utils.mod dotp_utils.mod kinds.mod \
                error_handling.mod timer.mod mp_interface.mod \
                system.mod parac.mod kpts.mod kpnt.mod spin.mod \
                elct.mod pslo.mod mm_extrap.mod legendre_p_utils.mod \
                ovlap_utils.mod rotate_utils.mod zeroing_utils.mod

mdfile_utils.mod.f90:$(SRCDIR)/mdfile_utils.mod.F90
mdfile_utils.mod.o:mdfile_utils.mod.f90 andp.mod andr.mod anneal_utils.mod \
                atwf.mod calc_alm_utils.mod cnst.mod comvel_utils.mod \
                comvelmod.mod coor.mod copot_utils.mod cotr.mod \
                ddipo_utils.mod detdof_utils.mod dispp_utils.mod \
                do_perturbation_p_utils.mod dynit_utils.mod \
                ekinpp_utils.mod elct.mod ener.mod error_handling.mod \
                extrap_utils.mod fileopen_utils.mod fileopenmod.mod \
                finalp_utils.mod fint.mod forcep_utils.mod \
                forces_diag_utils.mod geofile_utils.mod gsize_utils.mod \
                initrun_driver.mod initrun_utils.mod ions.mod \
                kinds.mod kpts.mod linres.mod localize_utils.mod \
                lr_tddft_utils.mod machine.mod mm_extrap.mod \
                moverho_utils.mod mp_interface.mod nlcc.mod \
                norm.mod nose.mod noseng_utils.mod nosepa_utils.mod \
                noseup_utils.mod nospinit_utils.mod parac.mod \
                phfac_utils.mod poin.mod posupi_utils.mod printave_utils.mod \
                printfor_utils.mod printp_utils.mod proppt_utils.mod \
                puttau_utils.mod rattle_utils.mod resetac_utils.mod \
                rhopri_utils.mod rinvel_utils.mod ropt.mod \
                rotvel_utils.mod rscvp_utils.mod sample_utils.mod \
                setirec_utils.mod shake_utils.mod soft.mod \
                spin.mod store_types.mod system.mod testex_utils.mod \
                teststore_utils.mod totstr_utils.mod vdwcmod.mod \
                velupi_utils.mod wrener_utils.mod wv30_utils.mod \
                zeroing_utils.mod

mdmain_utils.mod.f90:$(SRCDIR)/mdmain_utils.mod.F90
mdmain_utils.mod.o:mdmain_utils.mod.f90 anneal_utils.mod bicanonicalCpmd.mod \
                box_boundary_utils.mod bsym.mod bsympnt.mod \
                cnst.mod cnst_dyn.mod comvel_utils.mod comvelmod.mod \
                coor.mod copot_utils.mod cotr.mod csize_utils.mod \
                ddipo_utils.mod deort_utils.mod detdof_utils.mod \
                dispp_utils.mod dynit_utils.mod ekinpp_utils.mod \
                elct.mod ener.mod epr_efg_utils.mod error_handling.mod \
                fileopen_utils.mod fileopenmod.mod finalp_utils.mod \
                forcedr_driver.mod forcedr_utils.mod forcep_utils.mod \
                freqs_utils.mod geofile_utils.mod geq0mod.mod \
                gle_utils.mod glemod.mod gsize_utils.mod hfxmod.mod \
                hubbardu.mod initrun_driver.mod initrun_utils.mod \
                ions.mod isos.mod kinds.mod kpts.mod localize_utils.mod \
                lsforce_utils.mod machine.mod meta_colvar_inp_utils.mod \
                meta_colvar_utils.mod meta_exl_mult_utils.mod \
                meta_exlagr_methods.mod meta_exlagr_utils.mod \
                meta_multiple_walkers_utils.mod mfep.mod mp_interface.mod \
                mw.mod nlcc.mod norm.mod nose.mod noseinit_utils.mod \
                noseng_utils.mod nosepa_utils.mod noseup_utils.mod \
                nospinit_utils.mod ortho_utils.mod parac.mod \
                phfac_utils.mod pimd.mod posupa_utils.mod posupi_utils.mod \
                printave_utils.mod printp_utils.mod proja_utils.mod \
                pslo.mod puttau_utils.mod quenbo_utils.mod \
                rattle_utils.mod readsr_utils.mod rekine_utils.mod \
                resetac_utils.mod rhopri_utils.mod rinvel_utils.mod \
                ropt.mod rortv_utils.mod rotvel_utils.mod rscve_utils.mod \
                rscvp_utils.mod sample_utils.mod setbsstate_utils.mod \
                setirec_utils.mod shake_utils.mod soft.mod \
                spin.mod store_types.mod system.mod testex_utils.mod \
                teststore_utils.mod timer.mod totstr_utils.mod \
                tst2min_utils.mod utils.mod vdwcmod.mod velupa_utils.mod \
                velupi_utils.mod wann.mod wannier_print_utils.mod \
                wc_dos_utils.mod wrener_utils.mod wrgeo_utils.mod \
                wv30_utils.mod zeroing_utils.mod

mdpt_utils.mod.f90:$(SRCDIR)/mdpt_utils.mod.F90
mdpt_utils.mod.o:mdpt_utils.mod.f90 atwf.mod bsym.mod cl_init_utils.mod \
                clas.mod ddip.mod elct.mod error_handling.mod \
                fnlalloc_utils.mod fusion_utils.mod jrotation_utils.mod \
                kinds.mod kpts.mod linres.mod lr_in_utils.mod \
                md_driver.mod mdclas_utils.mod mdfile_utils.mod \
                mdmain_utils.mod mdshop_bo_utils.mod mdshop_cp_utils.mod \
                mm_dim_utils.mod mm_dimmod.mod mm_mddiag_utils.mod \
                mm_mdmain_utils.mod mm_mdshop_bo_utils.mod \
                mm_mdshop_cp_utils.mod mw.mod nabdy_md.mod \
                parac.mod prmem_utils.mod pslo.mod shop.mod \
                shop_rest_2.mod system.mod timer.mod utils.mod \
                vdwcmod.mod zeroing_utils.mod

mdshop_bo_utils.mod.f90:$(SRCDIR)/mdshop_bo_utils.mod.F90
mdshop_bo_utils.mod.o:mdshop_bo_utils.mod.f90 andp.mod andr.mod \
                anneal_utils.mod atwf.mod calc_alm_utils.mod \
                cnst.mod cnst_dyn.mod comvel_utils.mod comvelmod.mod \
                coor.mod copot_utils.mod cotr.mod deort_utils.mod \
                detdof_utils.mod dispp_utils.mod dynit_utils.mod \
                ekinpp_utils.mod elct.mod ener.mod error_handling.mod \
                extrap_utils.mod fileopen_utils.mod fileopenmod.mod \
                finalp_utils.mod fint.mod forcep_utils.mod \
                forces_diag_utils.mod geofile_utils.mod gsize_utils.mod \
                hubbardu.mod initrun_driver.mod ions.mod kinds.mod \
                kpts.mod linres.mod lr_tddft_utils.mod machine.mod \
                meta_colvar_inp_utils.mod meta_colvar_utils.mod \
                meta_exl_mult_utils.mod meta_exlagr_methods.mod \
                meta_exlagr_utils.mod mm_extrap.mod moverho_utils.mod \
                mp_interface.mod nlcc.mod norm.mod nose.mod \
                noseng_utils.mod nosepa_utils.mod noseup_utils.mod \
                nospinit_utils.mod parac.mod phfac_utils.mod \
                poin.mod posupi_utils.mod printave_utils.mod \
                printp_utils.mod proppt_utils.mod pslo.mod \
                puttau_utils.mod quenbo_utils.mod rattle_utils.mod \
                resetac_utils.mod response_pmod.mod rhoofr_utils.mod \
                rhopri_utils.mod rinitwf_utils.mod rinvel_utils.mod \
                rk4ov_utils.mod rnlsm_utils.mod ropt.mod rotvel_utils.mod \
                rscvp_utils.mod sample_utils.mod setirec_utils.mod \
                shake_utils.mod shop.mod shop_adds_utils.mod \
                shop_rest.mod shop_rest_2.mod soft.mod spin.mod \
                store_types.mod system.mod testex_utils.mod \
                teststore_utils.mod totstr_utils.mod velupi_utils.mod \
                wrener_utils.mod wrgeo_utils.mod wv30_utils.mod \
                zeroing_utils.mod

mdshop_cp_utils.mod.f90:$(SRCDIR)/mdshop_cp_utils.mod.F90
mdshop_cp_utils.mod.o:mdshop_cp_utils.mod.f90 anneal_utils.mod \
                cnst.mod comvel_utils.mod comvelmod.mod coor.mod \
                copot_utils.mod cotr.mod csize_utils.mod ddipo_utils.mod \
                deort_utils.mod detdof_utils.mod dispp_utils.mod \
                dynit_utils.mod ekinpp_utils.mod elct.mod ener.mod \
                error_handling.mod fileopen_utils.mod fileopenmod.mod \
                finalp_utils.mod forcedr_driver.mod forcedr_utils.mod \
                forcep_utils.mod geofile_utils.mod geq0mod.mod \
                gsize_utils.mod hubbardu.mod initrun_driver.mod \
                initrun_utils.mod kinds.mod machine.mod mp_interface.mod \
                nlcc.mod norm.mod nose.mod noseinit_utils.mod \
                noseng_utils.mod nosepa_utils.mod noseup_utils.mod \
                nospinit_utils.mod ortho_utils.mod parac.mod \
                phfac_utils.mod posupa_utils.mod posupi_utils.mod \
                printave_utils.mod printp_utils.mod proja_utils.mod \
                pslo.mod puttau_utils.mod quenbo_utils.mod \
                rattle_utils.mod rekine_utils.mod resetac_utils.mod \
                rhopri_utils.mod rinvel_utils.mod rk4ov_utils.mod \
                ropt.mod rortv_utils.mod rotvel_utils.mod rscve_utils.mod \
                rscvp_utils.mod sample_utils.mod setirec_utils.mod \
                shake_utils.mod shop.mod shop_adds_utils.mod \
                shop_rest.mod soft.mod spin.mod store_types.mod \
                system.mod testex_utils.mod teststore_utils.mod \
                totstr_utils.mod utils.mod velupa_utils.mod \
                velupi_utils.mod wann.mod wannier_print_utils.mod \
                wc_dos_utils.mod wrener_utils.mod wrgeo_utils.mod \
                wv30_utils.mod zeroing_utils.mod

merge.mod.f90:  $(SRCDIR)/merge.mod.F90
merge.mod.o:    merge.mod.f90 kinds.mod

merge_utils.mod.f90:$(SRCDIR)/merge_utils.mod.F90
merge_utils.mod.o:merge_utils.mod.f90 coor.mod cppt.mod elct.mod \
                fft_maxfft.mod filnmod.mod gsortho_utils.mod \
                kinds.mod mergemod.mod mp_interface.mod parac.mod \
                rhoofr_utils.mod rv30_utils.mod setirec_utils.mod \
                spin.mod store_types.mod system.mod

meta_cell_utils.mod.f90:$(SRCDIR)/meta_cell_utils.mod.F90
meta_cell_utils.mod.o:meta_cell_utils.mod.f90 cnst.mod cnst_dyn.mod \
                ener.mod error_handling.mod fileopen_utils.mod \
                fileopenmod.mod ions.mod kinds.mod latgen_utils.mod \
                meta_cv_utils.mod meta_hpot_utils.mod metr.mod \
                mp_interface.mod odiis_utils.mod parac.mod \
                prcp.mod readsr_utils.mod rmas.mod ropt.mod \
                soft.mod system.mod timer.mod tpar.mod zeroing_utils.mod

meta_colvar_inp_utils.mod.f90:$(SRCDIR)/meta_colvar_inp_utils.mod.F90
meta_colvar_inp_utils.mod.o:meta_colvar_inp_utils.mod.f90 chain_dr_utils.mod \
                cnst.mod cnst_dyn.mod cnstfc_utils.mod coninp_utils.mod \
                constr_utils.mod cotr.mod dum2_utils.mod error_handling.mod \
                fillc_utils.mod ions.mod kinds.mod meta_cv_qmmm_utils.mod \
                meta_cv_utils.mod meta_ex_mul_util_utils.mod \
                meta_exlagr_utils.mod meta_multiple_walkers_utils.mod \
                mm_dim_utils.mod mm_dimmod.mod mm_input.mod \
                mw.mod odiis_utils.mod parac.mod prcp.mod puttau_utils.mod \
                readsr_utils.mod system.mod timer.mod tpar.mod \
                zeroing_utils.mod

meta_colvar_utils.mod.f90:$(SRCDIR)/meta_colvar_utils.mod.F90
meta_colvar_utils.mod.o:meta_colvar_utils.mod.f90 cnst.mod \
                cnst_dyn.mod cotr.mod error_handling.mod fileopen_utils.mod \
                fileopenmod.mod ions.mod kinds.mod meta_colvar_inp_utils.mod \
                meta_colvar_util_utils.mod meta_exlagr_utils.mod \
                meta_hpot_utils.mod meta_multiple_walkers_utils.mod \
                mp_interface.mod mw.mod parac.mod pimd.mod \
                puttau_utils.mod ropt.mod soft.mod store_types.mod \
                strs.mod system.mod timer.mod tpar.mod zeroing_utils.mod

meta_colvar_util_utils.mod.f90:$(SRCDIR)/meta_colvar_util_utils.mod.F90
meta_colvar_util_utils.mod.o:meta_colvar_util_utils.mod.f90 \
                cnst.mod cnst_dyn.mod cotr.mod ener.mod error_handling.mod \
                fileopen_utils.mod fileopenmod.mod ions.mod \
                kinds.mod nose.mod parac.mod puttau_utils.mod \
                readsr_utils.mod ropt.mod system.mod tpar.mod \
                utils.mod zeroing_utils.mod

meta_cv_qmmm_utils.mod.f90:$(SRCDIR)/meta_cv_qmmm_utils.mod.F90
meta_cv_qmmm_utils.mod.o:meta_cv_qmmm_utils.mod.f90 cotr.mod \
                error_handling.mod fileopen_utils.mod fileopenmod.mod \
                fillc_utils.mod ions.mod isos.mod kinds.mod \
                meta_cv_utils.mod mm_dim_utils.mod mm_dimmod.mod \
                mm_input.mod parac.mod pbc_utils.mod system.mod \
                zeroing_utils.mod

meta_cv_utils.mod.f90:$(SRCDIR)/meta_cv_utils.mod.F90
meta_cv_utils.mod.o:meta_cv_utils.mod.f90 adat.mod cnst_dyn.mod \
                constr_utils.mod cotr.mod error_handling.mod \
                fileopen_utils.mod fileopenmod.mod fillc_utils.mod \
                ions.mod isos.mod kinds.mod latgen_utils.mod \
                meta_colvar_util_utils.mod metr.mod parac.mod \
                pbc_utils.mod system.mod zeroing_utils.mod

meta_dyn_def_utils.mod.f90:$(SRCDIR)/meta_dyn_def_utils.mod.F90
meta_dyn_def_utils.mod.o:meta_dyn_def_utils.mod.f90 cnst_dyn.mod \
                kinds.mod mw.mod system.mod

meta_exlagr_methods.mod.f90:$(SRCDIR)/meta_exlagr_methods.mod.F90
meta_exlagr_methods.mod.o:meta_exlagr_methods.mod.f90 cnst.mod \
                cnst_dyn.mod cotr.mod error_handling.mod fileopen_utils.mod \
                fileopenmod.mod ions.mod kinds.mod machine.mod \
                meta_colvar_inp_utils.mod meta_colvar_util_utils.mod \
                meta_exlagr_utils.mod meta_hpot_utils.mod meta_multiple_walkers_utils.mod \
                mm_dim_utils.mod mm_dimmod.mod mp_interface.mod \
                mw.mod parac.mod pimd.mod prng_utils.mod puttau_utils.mod \
                quenbo_utils.mod readsr_utils.mod ropt.mod \
                soft.mod strs.mod system.mod timer.mod tpar.mod \
                zeroing_utils.mod

meta_exlagr_utils.mod.f90:$(SRCDIR)/meta_exlagr_utils.mod.F90
meta_exlagr_utils.mod.o:meta_exlagr_utils.mod.f90 cnst.mod \
                cnst_dyn.mod ener.mod error_handling.mod fileopen_utils.mod \
                fileopenmod.mod ions.mod kinds.mod latgen_utils.mod \
                metr.mod mw.mod nose.mod parac.mod prng_utils.mod \
                readsr_utils.mod ropt.mod strs.mod system.mod \
                zeroing_utils.mod

meta_exl_mult_utils.mod.f90:$(SRCDIR)/meta_exl_mult_utils.mod.F90
meta_exl_mult_utils.mod.o:meta_exl_mult_utils.mod.f90 cnst.mod \
                cnst_dyn.mod cotr.mod error_handling.mod fileopen_utils.mod \
                fileopenmod.mod ions.mod kinds.mod machine.mod \
                meta_colvar_inp_utils.mod meta_ex_mul_util_utils.mod \
                meta_exlagr_utils.mod meta_hpot_utils.mod mp_interface.mod \
                parac.mod puttau_utils.mod readsr_utils.mod \
                ropt.mod soft.mod system.mod timer.mod tpar.mod \
                zeroing_utils.mod

meta_ex_mul_util_utils.mod.f90:$(SRCDIR)/meta_ex_mul_util_utils.mod.F90
meta_ex_mul_util_utils.mod.o:meta_ex_mul_util_utils.mod.f90 \
                cnst.mod cnst_dyn.mod ener.mod error_handling.mod \
                fileopen_utils.mod fileopenmod.mod kinds.mod \
                nose.mod parac.mod readsr_utils.mod ropt.mod \
                zeroing_utils.mod

metafun_utils.mod.f90:$(SRCDIR)/metafun_utils.mod.F90
metafun_utils.mod.o:metafun_utils.mod.f90 cnst.mod error_handling.mod \
                func.mod functionals_utils.mod kinds.mod lsd_func_utils.mod

meta_hpot_utils.mod.f90:$(SRCDIR)/meta_hpot_utils.mod.F90
meta_hpot_utils.mod.o:meta_hpot_utils.mod.f90 cnst.mod cnst_dyn.mod \
                error_handling.mod kinds.mod

meta_localizespin_utils.mod.f90:$(SRCDIR)/meta_localizespin_utils.mod.F90
meta_localizespin_utils.mod.o:meta_localizespin_utils.mod.f90 \
                cell.mod cnst_dyn.mod coor.mod fft_maxfft.mod \
                ions.mod kinds.mod latgen_utils.mod mp_interface.mod \
                parac.mod prcp.mod spin.mod system.mod

meta_multiple_walkers_utils.mod.f90:$(SRCDIR)/meta_multiple_walkers_utils.mod.F90
meta_multiple_walkers_utils.mod.o:meta_multiple_walkers_utils.mod.f90 \
                cnst.mod cnst_dyn.mod ener.mod error_handling.mod \
                fileopen_utils.mod fileopenmod.mod filnmod.mod \
                kinds.mod mp_interface.mod mw.mod nose.mod \
                parac.mod pimd.mod prng_utils.mod readsr_utils.mod \
                ropt.mod set_cp_grp_utils.mod system.mod timer.mod \
                wann.mod zeroing_utils.mod

metr.mod.f90:   $(SRCDIR)/metr.mod.F90
metr.mod.o:     metr.mod.f90 kinds.mod

mfep.mod.f90:   $(SRCDIR)/mfep.mod.F90
mfep.mod.o:     mfep.mod.f90 kinds.mod

min_heap.mod.f90:$(SRCDIR)/min_heap.mod.F90
min_heap.mod.o: min_heap.mod.f90 error_handling.mod

mixing_g_utils.mod.f90:$(SRCDIR)/mixing_g_utils.mod.F90
mixing_g_utils.mod.o:mixing_g_utils.mod.f90 anderson_utils.mod \
                andr.mod broy.mod broyden_utils.mod cppt.mod \
                error_handling.mod fftmain_utils.mod fftnew_utils.mod \
                geq0mod.mod kinds.mod parac.mod rhodiis_utils.mod \
                spin.mod system.mod zeroing_utils.mod

mixing_r_utils.mod.f90:$(SRCDIR)/mixing_r_utils.mod.F90
mixing_r_utils.mod.o:mixing_r_utils.mod.f90 anderson_utils.mod \
                andr.mod error_handling.mod kinds.mod rhodiis_utils.mod \
                spin.mod system.mod zeroing_utils.mod

mltfft_utils.mod.f90:$(SRCDIR)/mltfft_utils.mod.F90
mltfft_utils.mod.o:mltfft_utils.mod.f90 cublas_types.mod cublas_utils.mod \
                cuda_types.mod cufft_interfaces.mod cufft_types.mod \
                cufft_utils.mod cuuser_utils.mod error_handling.mod \
                gfft_utils.mod kinds.mod utils.mod timer.mod

mm_cpmd_add_MM_forces_f77_utils.mod.f90:$(SRCDIR)/mm_cpmd_add_MM_forces_f77_utils.mod.F90
mm_cpmd_add_MM_forces_f77_utils.mod.o:mm_cpmd_add_MM_forces_f77_utils.mod.f90 \
                ions.mod kinds.mod system.mod

mm_cpmd_esp_charges_f77_utils.mod.f90:$(SRCDIR)/mm_cpmd_esp_charges_f77_utils.mod.F90
mm_cpmd_esp_charges_f77_utils.mod.o:mm_cpmd_esp_charges_f77_utils.mod.f90 \
                cppt.mod elct.mod error_handling.mod fft_maxfft.mod \
                geq0mod.mod ions.mod isos.mod kinds.mod pslo.mod \
                system.mod zeroing_utils.mod

mm_cpmd_ext_pot_f77_utils.mod.f90:$(SRCDIR)/mm_cpmd_ext_pot_f77_utils.mod.F90
mm_cpmd_ext_pot_f77_utils.mod.o:mm_cpmd_ext_pot_f77_utils.mod.f90 \
                error_handling.mod kinds.mod parac.mod zeroing_utils.mod

mm_dim.mod.f90: $(SRCDIR)/mm_dim.mod.F90
mm_dim.mod.o:   mm_dim.mod.f90 kinds.mod

mm_dim_utils.mod.f90:$(SRCDIR)/mm_dim_utils.mod.F90
mm_dim_utils.mod.o:mm_dim_utils.mod.f90 error_handling.mod \
                ions.mod mm_dimmod.mod mm_input.mod parac.mod \
                system.mod timer.mod

mm_extrap.mod.f90:$(SRCDIR)/mm_extrap.mod.F90
mm_extrap.mod.o:mm_extrap.mod.f90 kinds.mod

mm_forcematch_utils.mod.f90:$(SRCDIR)/mm_forcematch_utils.mod.F90
mm_forcematch_utils.mod.o:mm_forcematch_utils.mod.f90 andp.mod \
                coor.mod copot_utils.mod efld.mod elct.mod \
                ener.mod error_handling.mod fileopen_utils.mod \
                fileopenmod.mod forcematch.mod forcematch_kfit_utils.mod \
                forcematch_qfit_utils.mod forcematch_utils.mod \
                forcep_utils.mod initrun_driver.mod ions.mod \
                kinds.mod kpts.mod machine.mod md_driver.mod \
                mm_dim_utils.mod mm_dimmod.mod mm_forces_diag_utils.mod \
                mm_input.mod mm_mddiag_utils.mod mm_parallel.mod \
                mp_interface.mod nlcc.mod parac.mod phfac_utils.mod \
                poin.mod prmem_utils.mod pslo.mod rinitwf_driver.mod \
                ropt.mod setirec_utils.mod spin.mod store_types.mod \
                system.mod wrener_utils.mod wrgeo_utils.mod \
                zeroing_utils.mod

mm_forces_diag_utils.mod.f90:$(SRCDIR)/mm_forces_diag_utils.mod.F90
mm_forces_diag_utils.mod.o:mm_forces_diag_utils.mod.f90 ehpsi_utils.mod \
                elct.mod ener.mod error_handling.mod kinds.mod \
                linres.mod lr_tddft_utils.mod machine.mod mm_extrap.mod \
                mm_input.mod mp_interface.mod norm.mod parac.mod \
                ropt.mod store_types.mod system.mod timer.mod \
                updrho_utils.mod updwf_utils.mod wrener_utils.mod \
                wv30_utils.mod zeroing_utils.mod

mm_forces_prop_utils.mod.f90:$(SRCDIR)/mm_forces_prop_utils.mod.F90
mm_forces_prop_utils.mod.o:mm_forces_prop_utils.mod.f90 ehpsi_utils.mod \
                ehrenfest_utils.mod elct.mod ener.mod error_handling.mod \
                geq0mod.mod k_updwf_utils.mod kinds.mod kpts.mod \
                linres.mod lr_tddft_utils.mod machine.mod mm_dim_utils.mod \
                mm_dimmod.mod mm_input.mod mp_interface.mod \
                parac.mod ropt.mod system.mod timer.mod updrho_utils.mod \
                updwf_utils.mod utils.mod zeroing_utils.mod

mm_init_utils.mod.f90:$(SRCDIR)/mm_init_utils.mod.F90
mm_init_utils.mod.o:mm_init_utils.mod.f90 adat.mod cell.mod \
                coor.mod cotr.mod error_handling.mod fileopen_utils.mod \
                fileopenmod.mod ions.mod isos.mod kinds.mod \
                mm_dimmod.mod mm_input.mod mm_ion_dens.mod \
                mm_parallel.mod mp_interface.mod parac.mod \
                rmas.mod store_types.mod system.mod zeroing_utils.mod

mm_input.mod.f90:$(SRCDIR)/mm_input.mod.F90
mm_input.mod.o: mm_input.mod.f90 kinds.mod

mm_ion_dens.mod.f90:$(SRCDIR)/mm_ion_dens.mod.F90
mm_ion_dens.mod.o:mm_ion_dens.mod.f90 kinds.mod system.mod

mm_mddiag_utils.mod.f90:$(SRCDIR)/mm_mddiag_utils.mod.F90
mm_mddiag_utils.mod.o:mm_mddiag_utils.mod.f90 andp.mod andr.mod \
                anneal_utils.mod atwf.mod calc_alm_utils.mod \
                cnst.mod cnst_dyn.mod comvel_utils.mod comvelmod.mod \
                coor.mod copot_utils.mod cotr.mod detdof_utils.mod \
                dispp_utils.mod dynit_utils.mod efld.mod ekinpp_utils.mod \
                elct.mod ener.mod error_handling.mod extrap_utils.mod \
                fileopen_utils.mod fileopenmod.mod finalp_utils.mod \
                fint.mod forcep_utils.mod forces_diag_utils.mod \
                geofile_utils.mod gsize_utils.mod hubbardu.mod \
                initrun_driver.mod initrun_utils.mod ions.mod \
                kinds.mod kpts.mod linres.mod localize_utils.mod \
                lr_tddft_utils.mod machine.mod md_driver.mod \
                mddiag_interaction_p_utils.mod meta_colvar_inp_utils.mod \
                meta_colvar_utils.mod meta_exl_mult_utils.mod \
                meta_exlagr_methods.mod mm_dim_utils.mod mm_dimmod.mod \
                mm_extrap.mod mm_forces_diag_utils.mod mm_forces_prop_utils.mod \
                mm_input.mod mm_mdmain_utils.mod mm_parallel.mod \
                moverho_utils.mod mp_interface.mod nlcc.mod \
                norm.mod nose.mod noseng_utils.mod nosepa_utils.mod \
                noseup_utils.mod nospinit_utils.mod parac.mod \
                phfac_utils.mod poin.mod posupi_utils.mod printave_utils.mod \
                printp_utils.mod prmem_utils.mod proppt_utils.mod \
                pslo.mod puttau_utils.mod rattle_utils.mod \
                resetac_utils.mod response_pmod.mod rhoofr_c_utils.mod \
                rhoofr_utils.mod rhopri_utils.mod rinitwf_utils.mod \
                rinvel_utils.mod rmas.mod rnlsm_utils.mod ropt.mod \
                rotvel_utils.mod rscvp_utils.mod sample_utils.mod \
                setirec_utils.mod sh_tddft_utils.mod shake_utils.mod \
                soc.mod soc_types.mod soft.mod spin.mod store_types.mod \
                system.mod td_input.mod td_utils.mod testex_utils.mod \
                teststore_utils.mod timer.mod totstr_utils.mod \
                tpar.mod tst2min_utils.mod vdwcmod.mod velupi_utils.mod \
                wrener_utils.mod wv30_utils.mod zeroing_utils.mod

mm_mdmain_utils.mod.f90:$(SRCDIR)/mm_mdmain_utils.mod.F90
mm_mdmain_utils.mod.o:mm_mdmain_utils.mod.f90 anneal_utils.mod \
                bsym.mod cnst.mod cnst_dyn.mod comvel_utils.mod \
                comvelmod.mod coor.mod copot_utils.mod cotr.mod \
                csize_utils.mod ddipo_utils.mod deort_utils.mod \
                detdof_utils.mod dispp_utils.mod dynit_utils.mod \
                efld.mod ekinpp_utils.mod elct.mod ener.mod \
                error_handling.mod fileopen_utils.mod fileopenmod.mod \
                finalp_utils.mod forcedr_driver.mod forcedr_utils.mod \
                forcep_utils.mod freqs_utils.mod geq0mod.mod \
                gsize_utils.mod hubbardu.mod initrun_driver.mod \
                initrun_utils.mod ions.mod kinds.mod kpts.mod \
                localize_utils.mod machine.mod meta_colvar_inp_utils.mod \
                meta_colvar_utils.mod meta_exl_mult_utils.mod \
                meta_exlagr_methods.mod meta_exlagr_utils.mod \
                mm_dim_utils.mod mm_dimmod.mod mm_input.mod \
                mm_parallel.mod mm_qmmm_forcedr_bs_utils.mod \
                mm_qmmm_forcedr_utils.mod mp_interface.mod \
                nlcc.mod norm.mod nose.mod noseinit_utils.mod \
                noseng_utils.mod nosepa_utils.mod noseup_utils.mod \
                nospinit_utils.mod ortho_utils.mod parac.mod \
                phfac_utils.mod posupa_utils.mod posupi_utils.mod \
                printave_utils.mod printp_utils.mod prmem_utils.mod \
                proja_utils.mod pslo.mod puttau_utils.mod quenbo_utils.mod \
                rattle_utils.mod rekine_utils.mod resetac_utils.mod \
                reshaper.mod rhopri_utils.mod rinvel_utils.mod \
                rmas.mod ropt.mod rortv_utils.mod rotvel_utils.mod \
                rscve_utils.mod rscvp_utils.mod sample_utils.mod \
                setbsstate_utils.mod setirec_utils.mod shake_utils.mod \
                shop_ekinqm.mod soft.mod spin.mod store_types.mod \
                system.mod testex_utils.mod teststore_utils.mod \
                timer.mod totstr_utils.mod tpar.mod tst2min_utils.mod \
                utils.mod vdwcmod.mod velupa_utils.mod velupi_utils.mod \
                wann.mod wannier_print_utils.mod wc_dos_utils.mod \
                wrener_utils.mod wrgeo_utils.mod wv30_utils.mod \
                zeroing_utils.mod

mm_mdshop_bo_utils.mod.f90:$(SRCDIR)/mm_mdshop_bo_utils.mod.F90
mm_mdshop_bo_utils.mod.o:mm_mdshop_bo_utils.mod.f90 andp.mod \
                andr.mod anneal_utils.mod atwf.mod calc_alm_utils.mod \
                cnst.mod cnst_dyn.mod comvel_utils.mod comvelmod.mod \
                coor.mod copot_utils.mod cotr.mod detdof_utils.mod \
                dispp_utils.mod dynit_utils.mod efld.mod ekinpp_utils.mod \
                elct.mod ener.mod error_handling.mod extrap_utils.mod \
                fileopen_utils.mod fileopenmod.mod finalp_utils.mod \
                fint.mod forcep_utils.mod forces_diag_utils.mod \
                geofile_utils.mod gsize_utils.mod hubbardu.mod \
                initrun_driver.mod ions.mod kinds.mod kpts.mod \
                linres.mod lr_tddft_utils.mod machine.mod mddiag_interaction_p_utils.mod \
                mdshop_bo_utils.mod meta_colvar_inp_utils.mod \
                meta_colvar_utils.mod meta_exl_mult_utils.mod \
                meta_exlagr_methods.mod mm_dim_utils.mod mm_dimmod.mod \
                mm_extrap.mod mm_forces_diag_utils.mod mm_input.mod \
                mm_mdmain_utils.mod mm_parallel.mod moverho_utils.mod \
                mp_interface.mod nlcc.mod norm.mod nose.mod \
                noseng_utils.mod nosepa_utils.mod noseup_utils.mod \
                nospinit_utils.mod parac.mod phfac_utils.mod \
                poin.mod posupi_utils.mod printave_utils.mod \
                printp_utils.mod prmem_utils.mod proppt_utils.mod \
                pslo.mod puttau_utils.mod rattle_utils.mod \
                resetac_utils.mod response_pmod.mod rhoofr_utils.mod \
                rhopri_utils.mod rinitwf_utils.mod rinvel_utils.mod \
                rk4ov_utils.mod rmas.mod rnlsm_utils.mod ropt.mod \
                rotvel_utils.mod rscvp_utils.mod sample_utils.mod \
                setirec_utils.mod shake_utils.mod shop.mod \
                shop_adds_utils.mod shop_rest.mod shop_rest_2.mod \
                soft.mod spin.mod store_types.mod system.mod \
                testex_utils.mod teststore_utils.mod totstr_utils.mod \
                tpar.mod tst2min_utils.mod velupi_utils.mod \
                wrener_utils.mod wrgeo_utils.mod wv30_utils.mod \
                zeroing_utils.mod

mm_mdshop_cp_utils.mod.f90:$(SRCDIR)/mm_mdshop_cp_utils.mod.F90
mm_mdshop_cp_utils.mod.o:mm_mdshop_cp_utils.mod.f90 anneal_utils.mod \
                cnst.mod cnst_dyn.mod comvel_utils.mod comvelmod.mod \
                coor.mod copot_utils.mod cotr.mod csize_utils.mod \
                ddipo_utils.mod deort_utils.mod detdof_utils.mod \
                dispp_utils.mod dynit_utils.mod efld.mod ekinpp_utils.mod \
                elct.mod ener.mod error_handling.mod fileopen_utils.mod \
                fileopenmod.mod finalp_utils.mod forcedr_driver.mod \
                forcedr_utils.mod forcep_utils.mod freqs_utils.mod \
                geq0mod.mod gsize_utils.mod hubbardu.mod initrun_driver.mod \
                initrun_utils.mod ions.mod kinds.mod kpts.mod \
                machine.mod meta_colvar_inp_utils.mod meta_colvar_utils.mod \
                meta_exl_mult_utils.mod meta_exlagr_methods.mod \
                meta_exlagr_utils.mod mm_dim_utils.mod mm_dimmod.mod \
                mm_input.mod mm_mdmain_utils.mod mm_parallel.mod \
                mm_qmmm_forcedr_utils.mod mp_interface.mod \
                nlcc.mod norm.mod nose.mod noseinit_utils.mod \
                noseng_utils.mod nosepa_utils.mod noseup_utils.mod \
                nospinit_utils.mod ortho_utils.mod parac.mod \
                phfac_utils.mod posupa_utils.mod posupi_utils.mod \
                printave_utils.mod printp_utils.mod proja_utils.mod \
                pslo.mod puttau_utils.mod quenbo_utils.mod \
                rattle_utils.mod rekine_utils.mod resetac_utils.mod \
                reshaper.mod rhopri_utils.mod rinvel_utils.mod \
                rk4ov_utils.mod rmas.mod ropt.mod rortv_utils.mod \
                rotvel_utils.mod rscve_utils.mod rscvp_utils.mod \
                sample_utils.mod setirec_utils.mod shake_utils.mod \
                shop.mod shop_adds_utils.mod shop_rest.mod \
                soft.mod spin.mod store_types.mod system.mod \
                testex_utils.mod teststore_utils.mod totstr_utils.mod \
                tpar.mod tst2min_utils.mod utils.mod velupa_utils.mod \
                velupi_utils.mod wann.mod wannier_print_utils.mod \
                wc_dos_utils.mod wrener_utils.mod wrgeo_utils.mod \
                wv30_utils.mod zeroing_utils.mod

mm_parallel.mod.f90:$(SRCDIR)/mm_parallel.mod.F90
mm_parallel.mod.o:mm_parallel.mod.f90

mm_qmmm_forcedr_bs_utils.mod.f90:$(SRCDIR)/mm_qmmm_forcedr_bs_utils.mod.F90
mm_qmmm_forcedr_bs_utils.mod.o:mm_qmmm_forcedr_bs_utils.mod.f90 \
                bsym.mod elct.mod ener.mod error_handling.mod \
                kinds.mod mm_dim_utils.mod mm_dimmod.mod mm_qmmm_forcedr_utils.mod \
                mp_interface.mod parac.mod setbsstate_utils.mod \
                system.mod zeroing_utils.mod

mm_qmmm_forcedr_utils.mod.f90:$(SRCDIR)/mm_qmmm_forcedr_utils.mod.F90
mm_qmmm_forcedr_utils.mod.o:mm_qmmm_forcedr_utils.mod.f90 bsym.mod \
                ener.mod error_handling.mod forcedr_driver.mod \
                kinds.mod machine.mod mm_dim_utils.mod mm_dimmod.mod \
                mm_input.mod mm_parallel.mod mp_interface.mod \
                norhoe_utils.mod parac.mod puttau_utils.mod \
                rhoofr_utils.mod rswfmod.mod system.mod timer.mod \
                utils.mod zeroing_utils.mod

mm_rho_forcedr_utils.mod.f90:$(SRCDIR)/mm_rho_forcedr_utils.mod.F90
mm_rho_forcedr_utils.mod.o:mm_rho_forcedr_utils.mod.f90 ener.mod \
                error_handling.mod kinds.mod machine.mod mm_input.mod \
                mm_parallel.mod mp_interface.mod parac.mod \
                puttau_utils.mod system.mod zeroing_utils.mod

molorb_utils.mod.f90:$(SRCDIR)/molorb_utils.mod.F90
molorb_utils.mod.o:molorb_utils.mod.f90 adat.mod cnst.mod empf.mod \
                empfor_utils.mod error_handling.mod ions.mod \
                kinds.mod mp_interface.mod ovlap_utils.mod \
                parac.mod pbc_utils.mod rotate_utils.mod spin.mod \
                system.mod wann.mod zeroing_utils.mod

mols.mod.f90:   $(SRCDIR)/mols.mod.F90
mols.mod.o:     mols.mod.f90

molstates_utils.mod.f90:$(SRCDIR)/molstates_utils.mod.F90
molstates_utils.mod.o:molstates_utils.mod.f90 cnst.mod ddip.mod \
                ddipo_utils.mod error_handling.mod forcep_utils.mod \
                hpsi_utils.mod ions.mod kinds.mod linres.mod \
                localize_utils.mod mm_input.mod mols.mod mp_interface.mod \
                ovlap_utils.mod parac.mod pbc_utils.mod poin.mod \
                rotate_utils.mod spin.mod system.mod utils.mod \
                wann.mod

molsym_utils.mod.f90:$(SRCDIR)/molsym_utils.mod.F90
molsym_utils.mod.o:molsym_utils.mod.f90 error_handling.mod \
                kinds.mod

moverho_utils.mod.f90:$(SRCDIR)/moverho_utils.mod.F90
moverho_utils.mod.o:moverho_utils.mod.f90 atrho_utils.mod atwf.mod \
                coor.mod cppt.mod error_handling.mod fft_maxfft.mod \
                fftmain_utils.mod fftnew_utils.mod fitpack_utils.mod \
                geq0mod.mod gvec.mod ions.mod kinds.mod parac.mod \
                prmem_utils.mod qspl.mod ropt.mod sfac.mod \
                spin.mod system.mod timer.mod zeroing_utils.mod

movi.mod.f90:   $(SRCDIR)/movi.mod.F90
movi.mod.o:     movi.mod.f90 system.mod

mp_interface.mod.f90:$(SRCDIR)/mp_interface.mod.F90
mp_interface.mod.o:mp_interface.mod.f90 error_handling.mod \
                kinds.mod machine.mod para_global.mod parac.mod \
                pstat.mod system.mod zeroing_utils.mod

mp_multiple_comm_init.mod.f90:$(SRCDIR)/mp_multiple_comm_init.mod.F90
mp_multiple_comm_init.mod.o:mp_multiple_comm_init.mod.f90 error_handling.mod \
                fileopen_utils.mod fileopenmod.mod kinds.mod \
                mp_interface.mod parac.mod pimd.mod readsr_utils.mod \
                set_cp_grp_utils.mod system.mod timer.mod zeroing_utils.mod

mtin_utils.mod.f90:$(SRCDIR)/mtin_utils.mod.F90
mtin_utils.mod.o:mtin_utils.mod.f90 cnst.mod cppt.mod error_handling.mod \
                fft_maxfft.mod fftmain_utils.mod geq0mod.mod \
                isos.mod kinds.mod parac.mod prmem_utils.mod \
                special_functions.mod system.mod timer.mod \
                zeroing_utils.mod

mulliken_utils.mod.f90:$(SRCDIR)/mulliken_utils.mod.F90
mulliken_utils.mod.o:mulliken_utils.mod.f90 atwf.mod augchg_utils.mod \
                dotp_utils.mod elct.mod error_handling.mod \
                ions.mod kinds.mod mp_interface.mod ovlap_utils.mod \
                parac.mod prop.mod pslo.mod setbasis_utils.mod \
                sfac.mod summat_utils.mod system.mod timer.mod \
                utils.mod zeroing_utils.mod

multtb_utils.mod.f90:$(SRCDIR)/multtb_utils.mod.F90
multtb_utils.mod.o:multtb_utils.mod.f90 error_handling.mod \
                kinds.mod parac.mod symm.mod

mw.mod.f90:     $(SRCDIR)/mw.mod.F90
mw.mod.o:       mw.mod.f90

my_para.f90:    $(SRCDIR)/my_para.F90
my_para.o:      my_para.f90 kinds.mod error_handling.mod timer.mod \
                mp_interface.mod mp_interface.mod kinds.mod \
                error_handling.mod timer.mod mp_interface.mod \
                utils.mod kinds.mod error_handling.mod timer.mod \
                mp_interface.mod kinds.mod error_handling.mod \
                timer.mod mp_interface.mod kinds.mod error_handling.mod \
                timer.mod machine.mod mp_interface.mod system.mod \
                parac.mod pstat.mod utils.mod kinds.mod error_handling.mod \
                timer.mod mp_interface.mod kinds.mod error_handling.mod \
                timer.mod mp_interface.mod parac.mod pstat.mod

nabdy_ampli.mod.f90:$(SRCDIR)/nabdy_ampli.mod.F90
nabdy_ampli.mod.o:nabdy_ampli.mod.f90 coor.mod cppt.mod error_handling.mod \
                fft_maxfft.mod fftmain_utils.mod geq0mod.mod \
                ions.mod kinds.mod mp_interface.mod nabdy_initialize.mod \
                nabdy_types.mod parac.mod prng_utils.mod rmas.mod \
                simulmod.mod system.mod zeroing_utils.mod

nabdy_forces.mod.f90:$(SRCDIR)/nabdy_forces.mod.F90
nabdy_forces.mod.o:nabdy_forces.mod.f90 coor.mod error_handling.mod \
                ions.mod kinds.mod nabdy_types.mod parac.mod \
                rmas.mod system.mod

nabdy_initialize.mod.f90:$(SRCDIR)/nabdy_initialize.mod.F90
nabdy_initialize.mod.o:nabdy_initialize.mod.f90 adat.mod cnst.mod \
                coor.mod error_handling.mod ions.mod kinds.mod \
                mm_dim_utils.mod mm_dimmod.mod mp_interface.mod \
                nabdy_types.mod parac.mod prng_utils.mod rmas.mod \
                system.mod zeroing_utils.mod

nabdy_md.mod.f90:$(SRCDIR)/nabdy_md.mod.F90
nabdy_md.mod.o: nabdy_md.mod.f90 andp.mod andr.mod anneal_utils.mod \
                atwf.mod box_boundary_utils.mod bs_forces_diag_utils.mod \
                bsym.mod calc_alm_utils.mod cdft_utils.mod \
                cdftmod.mod cnst.mod cnst_dyn.mod comvel_utils.mod \
                comvelmod.mod coor.mod copot_utils.mod cotr.mod \
                detdof_utils.mod dispp_utils.mod dynit_utils.mod \
                ekinpp_utils.mod elct.mod ener.mod error_handling.mod \
                extrap_utils.mod fileopen_utils.mod fileopenmod.mod \
                finalp_utils.mod fint.mod forcep_utils.mod \
                forces_diag_utils.mod geofile_utils.mod gle_utils.mod \
                glemod.mod gsize_utils.mod hfxmod.mod initrun_driver.mod \
                ions.mod isos.mod kinds.mod kpts.mod linres.mod \
                localize_utils.mod lr_tddft_utils.mod machine.mod \
                md_driver.mod mddiag_interaction_p_utils.mod \
                meta_multiple_walkers_utils.mod mfep.mod mm_dim_utils.mod \
                mm_dimmod.mod mm_extrap.mod moverho_utils.mod \
                mp_interface.mod nabdy_ampli.mod nabdy_forces.mod \
                nabdy_initialize.mod nabdy_types.mod nlcc.mod \
                norm.mod nose.mod noseng_utils.mod nosepa_utils.mod \
                noseup_utils.mod nospinit_utils.mod parac.mod \
                phfac_utils.mod pimd.mod poin.mod posupi_utils.mod \
                printave_utils.mod printp_utils.mod prng_utils.mod \
                proppt_utils.mod puttau_utils.mod rattle_utils.mod \
                readsr_utils.mod resetac_utils.mod response_pmod.mod \
                rhopri_utils.mod rinvel_utils.mod rmas.mod \
                ropt.mod rotvel_utils.mod rscvp_utils.mod sample_utils.mod \
                setbsstate_utils.mod setirec_utils.mod shake_utils.mod \
                soft.mod spin.mod store_types.mod system.mod \
                testex_utils.mod teststore_utils.mod timer.mod \
                totstr_utils.mod tpar.mod vdwcmod.mod velupi_utils.mod \
                wrener_utils.mod wrgeo_utils.mod wv30_utils.mod \
                zeroing_utils.mod

nabdy_types.mod.f90:$(SRCDIR)/nabdy_types.mod.F90
nabdy_types.mod.o:nabdy_types.mod.f90 kinds.mod

newcell_utils.mod.f90:$(SRCDIR)/newcell_utils.mod.F90
newcell_utils.mod.o:newcell_utils.mod.f90 initclust_utils.mod \
                nlccset_utils.mod rggen_utils.mod rinforce_utils.mod

newd_utils.mod.f90:$(SRCDIR)/newd_utils.mod.F90
newd_utils.mod.o:newd_utils.mod.f90 cppt.mod error_handling.mod \
                geq0mod.mod ions.mod kinds.mod mp_interface.mod \
                nlps.mod parac.mod pslo.mod qvan2_utils.mod \
                sfac.mod system.mod timer.mod zeroing_utils.mod \
                kinds.mod error_handling.mod timer.mod

nfunc_utils.mod.f90:$(SRCDIR)/nfunc_utils.mod.F90
nfunc_utils.mod.o:nfunc_utils.mod.f90 error_handling.mod kinds.mod \
                timer.mod error_handling.mod kinds.mod timer.mod \
                error_handling.mod kinds.mod timer.mod error_handling.mod \
                kinds.mod timer.mod error_handling.mod kinds.mod \
                timer.mod error_handling.mod kinds.mod timer.mod \
                error_handling.mod kinds.mod timer.mod error_handling.mod \
                kinds.mod timer.mod error_handling.mod kinds.mod \
                timer.mod error_handling.mod kinds.mod timer.mod \
                error_handling.mod kinds.mod timer.mod error_handling.mod \
                kinds.mod timer.mod

nlcc.mod.f90:   $(SRCDIR)/nlcc.mod.F90
nlcc.mod.o:     nlcc.mod.f90 kinds.mod system.mod

nlccset_utils.mod.f90:$(SRCDIR)/nlccset_utils.mod.F90
nlccset_utils.mod.o:nlccset_utils.mod.f90 bessm_utils.mod cnst.mod \
                cppt.mod error_handling.mod fitpack_utils.mod \
                ions.mod kinds.mod mp_interface.mod nlcc.mod \
                parac.mod pslo.mod qspl.mod radin_utils.mod \
                system.mod timer.mod vdbp.mod zeroing_utils.mod

nlccstr_utils.mod.f90:$(SRCDIR)/nlccstr_utils.mod.F90
nlccstr_utils.mod.o:nlccstr_utils.mod.f90 copot_utils.mod cppt.mod \
                dotp_utils.mod fftmain_utils.mod kinds.mod \
                nlcc.mod parac.mod sfac.mod strs.mod system.mod \
                timer.mod utils.mod

nlforce_utils.mod.f90:$(SRCDIR)/nlforce_utils.mod.F90
nlforce_utils.mod.o:nlforce_utils.mod.f90 cppt.mod cvan.mod \
                error_handling.mod ions.mod kinds.mod mp_interface.mod \
                nlps.mod parac.mod pslo.mod sfac.mod sgpp.mod \
                spin.mod system.mod tbxc.mod timer.mod zeroing_utils.mod

nlps.mod.f90:   $(SRCDIR)/nlps.mod.F90
nlps.mod.o:     nlps.mod.f90 kinds.mod system.mod

nl_res_utils.mod.f90:$(SRCDIR)/nl_res_utils.mod.F90
nl_res_utils.mod.o:nl_res_utils.mod.f90 cppt.mod elct.mod error_handling.mod \
                kinds.mod mp_interface.mod nlps.mod parac.mod \
                pslo.mod sfac.mod sgpp.mod system.mod timer.mod \
                zeroing_utils.mod

nlsl_utils.mod.f90:$(SRCDIR)/nlsl_utils.mod.F90
nlsl_utils.mod.o:nlsl_utils.mod.f90 cvan.mod error_handling.mod \
                ions.mod kinds.mod nlps.mod parac.mod pslo.mod \
                sfac.mod str2.mod strs.mod system.mod timer.mod

nlsm1_s_utils.mod.f90:$(SRCDIR)/nlsm1_s_utils.mod.F90
nlsm1_s_utils.mod.o:nlsm1_s_utils.mod.f90 geq0mod.mod ions.mod \
                kinds.mod kpnt.mod kpts.mod mp_interface.mod \
                nlps.mod parac.mod str2.mod system.mod timer.mod

nmr_chi_p_utils.mod.f90:$(SRCDIR)/nmr_chi_p_utils.mod.F90
nmr_chi_p_utils.mod.o:nmr_chi_p_utils.mod.f90 elct.mod error_handling.mod \
                kinds.mod mp_interface.mod nmr_position_p_utils.mod \
                nmr_util_p_utils.mod parac.mod response_pmod.mod \
                system.mod timer.mod zeroing_utils.mod

nmr_current_p_utils.mod.f90:$(SRCDIR)/nmr_current_p_utils.mod.F90
nmr_current_p_utils.mod.o:nmr_current_p_utils.mod.f90 coor.mod \
                cppt.mod dotp_utils.mod elct.mod error_handling.mod \
                forcep_utils.mod geq0mod.mod ions.mod kinds.mod \
                mp_interface.mod nmr_position_p_utils.mod nmr_util_p_utils.mod \
                parac.mod prmem_utils.mod response_pmod.mod \
                system.mod timer.mod zeroing_utils.mod

nmr_full_p_utils.mod.f90:$(SRCDIR)/nmr_full_p_utils.mod.F90
nmr_full_p_utils.mod.o:nmr_full_p_utils.mod.f90 error_handling.mod \
                kinds.mod nmr_position_p_utils.mod nmr_util_p_utils.mod \
                parac.mod response_pmod.mod system.mod timer.mod \
                zeroing_utils.mod kinds.mod error_handling.mod \
                timer.mod system.mod parac.mod soft.mod response_pmod.mod \
                nmr_util_p_utils.mod nmr_shift_p_utils.mod \
                nmr_chi_p_utils.mod nmr_util_p_utils.mod nmr_position_p_utils.mod \
                perturbation_p_utils.mod rwfopt_p_utils.mod \
                zeroing_utils.mod fft_maxfft.mod

nmr_para_p_utils.mod.f90:$(SRCDIR)/nmr_para_p_utils.mod.F90
nmr_para_p_utils.mod.o:nmr_para_p_utils.mod.f90 error_handling.mod \
                kinds.mod mp_interface.mod parac.mod response_pmod.mod \
                system.mod

nmr_position_p_utils.mod.f90:$(SRCDIR)/nmr_position_p_utils.mod.F90
nmr_position_p_utils.mod.o:nmr_position_p_utils.mod.f90 coor.mod \
                dotp_utils.mod error_handling.mod fft_maxfft.mod \
                geq0mod.mod gvec.mod ions.mod kinds.mod mp_interface.mod \
                nmr_util_p_utils.mod parac.mod reshaper.mod \
                response_pmod.mod system.mod timer.mod util_p_utils.mod \
                zeroing_utils.mod kinds.mod error_handling.mod \
                timer.mod system.mod parac.mod response_pmod.mod \
                nmr_position_p_utils.mod

nmr_p_utils.mod.f90:$(SRCDIR)/nmr_p_utils.mod.F90
nmr_p_utils.mod.o:nmr_p_utils.mod.f90 cnst.mod coor.mod ddip.mod \
                ddipo_utils.mod eicalc_utils.mod elct.mod error_handling.mod \
                fft_maxfft.mod isos.mod kinds.mod localize_utils.mod \
                machine.mod mp_interface.mod nmr_chi_p_utils.mod \
                nmr_current_p_utils.mod nmr_full_p_utils.mod \
                nmr_para_p_utils.mod nmr_position_p_utils.mod \
                nmr_shift_p_utils.mod nmr_util_p_utils.mod \
                parac.mod perturbation_p_utils.mod phfac_utils.mod \
                prmem_utils.mod prop.mod response_pmod.mod \
                restart_p_utils.mod ropt.mod rwfopt_p_utils.mod \
                setirec_utils.mod soft.mod store_types.mod \
                system.mod timer.mod utils.mod wann.mod wv30_utils.mod \
                zeroing_utils.mod

nmr_shift_p_utils.mod.f90:$(SRCDIR)/nmr_shift_p_utils.mod.F90
nmr_shift_p_utils.mod.o:nmr_shift_p_utils.mod.f90 adat.mod \
                cnst.mod cppt.mod elct.mod error_handling.mod \
                fft_maxfft.mod geq0mod.mod ions.mod kinds.mod \
                mp_interface.mod nmr_position_p_utils.mod nmr_util_p_utils.mod \
                parac.mod response_pmod.mod system.mod timer.mod \
                zeroing_utils.mod kinds.mod error_handling.mod \
                timer.mod parac.mod response_pmod.mod ions.mod \
                kinds.mod error_handling.mod timer.mod system.mod \
                parac.mod cnst.mod geq0mod.mod cppt.mod

nmr_util_p_utils.mod.f90:$(SRCDIR)/nmr_util_p_utils.mod.F90
nmr_util_p_utils.mod.o:nmr_util_p_utils.mod.f90 coor.mod error_handling.mod \
                fftmain_utils.mod fftutil_utils.mod ions.mod \
                kinds.mod machine.mod parac.mod response_pmod.mod \
                system.mod timer.mod zeroing_utils.mod kinds.mod \
                error_handling.mod timer.mod mp_interface.mod \
                system.mod parac.mod cppt.mod response_pmod.mod \
                fftutil_utils.mod fftmain_utils.mod fft_maxfft.mod \
                zeroing_utils.mod kinds.mod error_handling.mod \
                timer.mod mp_interface.mod system.mod parac.mod \
                cppt.mod response_pmod.mod fftutil_utils.mod \
                fftmain_utils.mod zeroing_utils.mod fft_maxfft.mod \
                kinds.mod error_handling.mod timer.mod system.mod \
                parac.mod coor.mod ions.mod response_pmod.mod \
                fft_maxfft.mod zeroing_utils.mod kinds.mod \
                error_handling.mod timer.mod system.mod parac.mod \
                cppt.mod fftutil_utils.mod fftmain_utils.mod \
                zeroing_utils.mod fft_maxfft.mod kinds.mod \
                error_handling.mod timer.mod system.mod parac.mod \
                cppt.mod kinds.mod error_handling.mod timer.mod \
                system.mod parac.mod kinds.mod error_handling.mod \
                timer.mod system.mod parac.mod cppt.mod fftutil_utils.mod \
                fftmain_utils.mod zeroing_utils.mod fft_maxfft.mod

nofo.mod.f90:   $(SRCDIR)/nofo.mod.F90
nofo.mod.o:     nofo.mod.f90 kinds.mod

noforce_utils.mod.f90:$(SRCDIR)/noforce_utils.mod.F90
noforce_utils.mod.o:noforce_utils.mod.f90 csmat_utils.mod dotp_utils.mod \
                elct.mod error_handling.mod fnonloc_utils.mod \
                geq0mod.mod gsize_utils.mod hnlmat_utils.mod \
                hubbardu.mod hubbardu_utils.mod kinds.mod mp_interface.mod \
                nlforce_utils.mod nlps.mod nlsl_utils.mod norm.mod \
                ovlap_utils.mod parac.mod pslo.mod puttau_utils.mod \
                reigs_utils.mod rgs_utils.mod rnlfl_utils.mod \
                rnlsm_utils.mod ropt.mod rotate_utils.mod rscpot_utils.mod \
                sfac.mod spin.mod summat_utils.mod symtrz_utils.mod \
                system.mod timer.mod utils.mod vpsi_utils.mod \
                zeroing_utils.mod

norhoe_utils.mod.f90:$(SRCDIR)/norhoe_utils.mod.F90
norhoe_utils.mod.o:norhoe_utils.mod.f90 csmat_utils.mod error_handling.mod \
                fft_maxfft.mod geq0mod.mod jacobi_utils.mod \
                kinds.mod mp_interface.mod nlps.mod noforce_utils.mod \
                parac.mod pslo.mod rhoofr_utils.mod rnlsm_utils.mod \
                rotate_utils.mod sfac.mod spin.mod system.mod \
                timer.mod utils.mod

norm.mod.f90:   $(SRCDIR)/norm.mod.F90
norm.mod.o:     norm.mod.f90 kinds.mod

nort.mod.f90:   $(SRCDIR)/nort.mod.F90
nort.mod.o:     nort.mod.f90 kinds.mod

nosalloc_utils.mod.f90:$(SRCDIR)/nosalloc_utils.mod.F90
nosalloc_utils.mod.o:nosalloc_utils.mod.f90 bsym.mod error_handling.mod \
                fileopen_utils.mod fileopenmod.mod ions.mod \
                mm_dim_utils.mod mm_dimmod.mod mm_input.mod \
                mp_interface.mod mw.mod nose.mod parac.mod \
                pimd.mod prmem_utils.mod system.mod timer.mod

noscinit_utils.mod.f90:$(SRCDIR)/noscinit_utils.mod.F90
noscinit_utils.mod.o:noscinit_utils.mod.f90 cnst.mod nose.mod \
                system.mod zeroing_utils.mod

noseinit_utils.mod.f90:$(SRCDIR)/noseinit_utils.mod.F90
noseinit_utils.mod.o:noseinit_utils.mod.f90 kinds.mod nose.mod \
                parac.mod system.mod zeroing_utils.mod

nose.mod.f90:   $(SRCDIR)/nose.mod.F90
nose.mod.o:     nose.mod.f90 kinds.mod system.mod

noseng_utils.mod.f90:$(SRCDIR)/noseng_utils.mod.F90
noseng_utils.mod.o:noseng_utils.mod.f90 bsym.mod cnst.mod fileopen_utils.mod \
                fileopenmod.mod ions.mod kinds.mod nose.mod \
                parac.mod pimd.mod rmas.mod system.mod zeroing_utils.mod

nosepa_utils.mod.f90:$(SRCDIR)/nosepa_utils.mod.F90
nosepa_utils.mod.o:nosepa_utils.mod.f90 cnst.mod cotr.mod elct.mod \
                error_handling.mod ions.mod isos.mod kinds.mod \
                mm_dimmod.mod mm_input.mod nose.mod parac.mod \
                pimd.mod prcp.mod system.mod zeroing_utils.mod

noseup_utils.mod.f90:$(SRCDIR)/noseup_utils.mod.F90
noseup_utils.mod.o:noseup_utils.mod.f90 bsym.mod enosmove_utils.mod \
                error_handling.mod kinds.mod mp_interface.mod \
                nose.mod parac.mod pimd.mod pnosmove_utils.mod \
                prcnosmove_utils.mod prcp.mod prpcmove_utils.mod \
                prpcnosmove_utils.mod prpnosmove_utils.mod \
                rekine_utils.mod rmas.mod system.mod timer.mod

nospinit_utils.mod.f90:$(SRCDIR)/nospinit_utils.mod.F90
nospinit_utils.mod.o:nospinit_utils.mod.f90 cnst.mod ions.mod \
                kinds.mod nose.mod parac.mod pimd.mod prng_utils.mod \
                system.mod zeroing_utils.mod

npt_md_utils.mod.f90:$(SRCDIR)/npt_md_utils.mod.F90
npt_md_utils.mod.o:npt_md_utils.mod.f90 andp.mod andr.mod anneal_utils.mod \
                cnst.mod comvel_utils.mod comvelmod.mod coor.mod \
                copot_utils.mod csize_utils.mod ddipo_utils.mod \
                deort_utils.mod detdof_utils.mod dispp_utils.mod \
                dynit_utils.mod ekinpp_utils.mod elct.mod ener.mod \
                error_handling.mod extrap_utils.mod fileopen_utils.mod \
                fileopenmod.mod finalp_utils.mod forcedr_driver.mod \
                forcedr_utils.mod forcep_utils.mod forces_diag_utils.mod \
                freqs_utils.mod geofile_utils.mod geq0mod.mod \
                gsize_utils.mod initrun_driver.mod initrun_utils.mod \
                kinds.mod kpts.mod localize_utils.mod machine.mod \
                metr.mod mm_extrap.mod mp_interface.mod newcell_utils.mod \
                nlcc.mod norm.mod noscinit_utils.mod nose.mod \
                noseinit_utils.mod noseng_utils.mod nosepa_utils.mod \
                noseup_utils.mod nospinit_utils.mod ortho_utils.mod \
                parac.mod phfac_utils.mod poin.mod posupa_utils.mod \
                posupi_utils.mod prcp.mod printave_utils.mod \
                printp_utils.mod proja_utils.mod pslo.mod puttau_utils.mod \
                quenbo_utils.mod rattle_utils.mod rekine_utils.mod \
                resetac_utils.mod rhopri_utils.mod rinvel_utils.mod \
                ropt.mod rortv_utils.mod rscve_utils.mod rscvp_utils.mod \
                setirec_utils.mod setsc_utils.mod shock.mod \
                soft.mod spin.mod store_types.mod system.mod \
                testex_utils.mod teststore_utils.mod totstr_utils.mod \
                utils.mod vdwcmod.mod velupa_utils.mod velupi_utils.mod \
                vepsup_utils.mod wrener_utils.mod wrgeo_utils.mod \
                wv30_utils.mod zeroing_utils.mod

nuclear_p_utils.mod.f90:$(SRCDIR)/nuclear_p_utils.mod.F90
nuclear_p_utils.mod.o:nuclear_p_utils.mod.f90 cnst.mod coor.mod \
                densrd_utils.mod elct.mod error_handling.mod \
                fft_maxfft.mod ions.mod kinds.mod machine.mod \
                mm_dim_utils.mod mm_dimmod.mod mp_interface.mod \
                parac.mod printp_utils.mod recpnew_utils.mod \
                response_pmod.mod rhoofr_p_utils.mod rhoofr_utils.mod \
                rhopri_utils.mod rinforce_nuc_utils.mod ropt.mod \
                rwfopt_nuc_utils.mod rwfopt_p_utils.mod sgpp.mod \
                soft.mod spin.mod system.mod testex_utils.mod \
                write_pp_utils.mod zeroing_utils.mod

numpw_utils.mod.f90:$(SRCDIR)/numpw_utils.mod.F90
numpw_utils.mod.o:numpw_utils.mod.f90 cell.mod cnst.mod error_handling.mod \
                fftchk_utils.mod fint.mod gvec.mod kinds.mod \
                kpts.mod parac.mod rggen_utils.mod sort_utils.mod \
                sphe.mod system.mod timer.mod

nvar.mod.f90:   $(SRCDIR)/nvar.mod.F90
nvar.mod.o:     nvar.mod.f90

nvtx_interfaces.mod.f90:$(SRCDIR)/nvtx_interfaces.mod.F90
nvtx_interfaces.mod.o:nvtx_interfaces.mod.f90

nvtx_utils.mod.f90:$(SRCDIR)/nvtx_utils.mod.F90
nvtx_utils.mod.o:nvtx_utils.mod.f90 nvtx_interfaces.mod

odiis_p_utils.mod.f90:$(SRCDIR)/odiis_p_utils.mod.F90
odiis_p_utils.mod.o:odiis_p_utils.mod.f90 dotp_utils.mod ener.mod \
                kinds.mod mp_interface.mod odiis_utils.mod \
                parac.mod pcgrad_p_utils.mod system.mod timer.mod \
                zeroing_utils.mod

odiis_utils.mod.f90:$(SRCDIR)/odiis_utils.mod.F90
odiis_utils.mod.o:odiis_utils.mod.f90 error_handling.mod kinds.mod \
                parac.mod system.mod timer.mod zeroing_utils.mod \
                kinds.mod error_handling.mod timer.mod kinds.mod \
                error_handling.mod timer.mod kinds.mod error_handling.mod \
                timer.mod kinds.mod error_handling.mod timer.mod \
                kinds.mod error_handling.mod timer.mod mp_interface.mod \
                system.mod parac.mod kinds.mod error_handling.mod \
                timer.mod kinds.mod error_handling.mod timer.mod \
                mp_interface.mod system.mod parac.mod ener.mod \
                elct.mod cp_grp_utils.mod cp_grp_utils.mod \
                odiis_utils.mod zeroing_utils.mod nvtx_utils.mod

ohfd_utils.mod.f90:$(SRCDIR)/ohfd_utils.mod.F90
ohfd_utils.mod.o:ohfd_utils.mod.f90 andp.mod canon_utils.mod \
                coor.mod copot_utils.mod dynit_utils.mod elct.mod \
                error_handling.mod fileopen_utils.mod fileopenmod.mod \
                forcedr_driver.mod forcep_utils.mod forces_diag_utils.mod \
                initrun_driver.mod kinds.mod machine.mod mp_interface.mod \
                nlcc.mod parac.mod phfac_utils.mod poin.mod \
                pslo.mod rhoofr_utils.mod rinitwf_utils.mod \
                rnlsm_utils.mod ropt.mod setirec_utils.mod \
                soft.mod spin.mod store_types.mod system.mod \
                wv30_utils.mod zeroing_utils.mod

ohlr_utils.mod.f90:$(SRCDIR)/ohlr_utils.mod.F90
ohlr_utils.mod.o:ohlr_utils.mod.f90 atwf.mod canon_utils.mod \
                coor.mod cppt.mod ddipo_utils.mod dotp_utils.mod \
                dynit_utils.mod elct.mod ener.mod error_handling.mod \
                fftmain_utils.mod fileopen_utils.mod fileopenmod.mod \
                fnlalloc_utils.mod forcedr_driver.mod forcedr_utils.mod \
                forcep_utils.mod geq0mod.mod initrun_driver.mod \
                initrun_utils.mod ions.mod isos.mod kinds.mod \
                ksdiag_utils.mod linres.mod localize_utils.mod \
                lr_in_utils.mod lr_xcpot_utils.mod machine.mod \
                mp_interface.mod nlcc.mod norm.mod opt_lr_utils.mod \
                ovlap_utils.mod parac.mod pbc_utils.mod poin.mod \
                puttau_utils.mod rho1ofr_utils.mod rhoofr_utils.mod \
                rnlsm_utils.mod ropt.mod setbasis_utils.mod \
                setirec_utils.mod soft.mod spin.mod store_types.mod \
                system.mod testex_utils.mod timer.mod updwf_utils.mod \
                utils.mod v1ofrho1_utils.mod vhk_utils.mod \
                vpsi_utils.mod wrener_utils.mod wv30_utils.mod \
                zeroing_utils.mod kinds.mod mp_interface.mod \
                error_handling.mod timer.mod system.mod parac.mod

opeigr_c_utils.mod.f90:$(SRCDIR)/opeigr_c_utils.mod.F90
opeigr_c_utils.mod.o:opeigr_c_utils.mod.f90 ddip.mod error_handling.mod \
                kinds.mod mp_interface.mod opeigr_utils.mod \
                parac.mod reshaper.mod spin.mod system.mod \
                timer.mod utils.mod zeroing_utils.mod

opeigr_p_utils.mod.f90:$(SRCDIR)/opeigr_p_utils.mod.F90
opeigr_p_utils.mod.o:opeigr_p_utils.mod.f90 ddip.mod error_handling.mod \
                kinds.mod mp_interface.mod opeigr_utils.mod \
                parac.mod reshaper.mod spin.mod system.mod \
                zeroing_utils.mod

opeigr_utils.mod.f90:$(SRCDIR)/opeigr_utils.mod.F90
opeigr_utils.mod.o:opeigr_utils.mod.f90 ddip.mod error_handling.mod \
                gvec.mod kinds.mod kpts.mod mp_interface.mod \
                numpw_utils.mod parac.mod sphe.mod spin.mod \
                system.mod timer.mod utils.mod zeroing_utils.mod

opt_lr_utils.mod.f90:$(SRCDIR)/opt_lr_utils.mod.F90
opt_lr_utils.mod.o:opt_lr_utils.mod.f90 cppt.mod elct.mod error_handling.mod \
                fft.mod fftmain_utils.mod fftnew_utils.mod \
                geq0mod.mod kinds.mod linres.mod lr_ortho_utils.mod \
                lr_upd_utils.mod machine.mod norm.mod parac.mod \
                ropt.mod spin.mod system.mod timer.mod tpar.mod \
                utils.mod vpsi_utils.mod zeroing_utils.mod

orbhard_utils.mod.f90:$(SRCDIR)/orbhard_utils.mod.F90
orbhard_utils.mod.o:orbhard_utils.mod.f90 atwf.mod ddip.mod \
                elct.mod elct2.mod error_handling.mod fnlalloc_utils.mod \
                inscan_utils.mod kinds.mod kpts.mod linres.mod \
                ohfd_utils.mod ohlr_utils.mod parac.mod system.mod \
                utils.mod wann.mod

orbrot_utils.mod.f90:$(SRCDIR)/orbrot_utils.mod.F90
orbrot_utils.mod.o:orbrot_utils.mod.f90 error_handling.mod \
                kinds.mod mp_interface.mod parac.mod timer.mod \
                utils.mod zeroing_utils.mod

ortho_utils.mod.f90:$(SRCDIR)/ortho_utils.mod.F90
ortho_utils.mod.o:ortho_utils.mod.f90 cp_cuda_types.mod cp_cuortho_types.mod \
                cp_cuortho_utils.mod cp_cuwfn_types.mod cp_cuwfn_utils.mod \
                cp_grp_utils.mod disortho_utils.mod elct.mod \
                error_handling.mod geq0mod.mod gsortho_utils.mod \
                jrotation_utils.mod kinds.mod kpts.mod lowdin_utils.mod \
                mp_interface.mod parac.mod part_1d.mod pslo.mod \
                rgs_utils.mod rgsvan_utils.mod sfac.mod spin.mod \
                system.mod td_input.mod timer.mod utils.mod

ovlap_utils.mod.f90:$(SRCDIR)/ovlap_utils.mod.F90
ovlap_utils.mod.o:ovlap_utils.mod.f90 cp_cuwfn_types.mod cp_grp_utils.mod \
                cublas_types.mod cublas_utils.mod cuda_types.mod \
                cuda_utils.mod error_handling.mod geq0mod.mod \
                kinds.mod nvtx_utils.mod parac.mod sizeof_kinds.mod \
                spin.mod system.mod timer.mod zeroing_utils.mod \
                system.mod cp_grp_utils.mod jrotation_utils.mod \
                mp_interface.mod kinds.mod error_handling.mod \
                timer.mod zeroing_utils.mod kinds.mod error_handling.mod \
                timer.mod system.mod parac.mod zeroing_utils.mod \
                kinds.mod error_handling.mod timer.mod geq0mod.mod \
                cp_grp_utils.mod cp_grp_utils.mod zeroing_utils.mod \
                kinds.mod error_handling.mod timer.mod system.mod \
                parac.mod spin.mod cp_grp_utils.mod cp_grp_utils.mod \
                zeroing_utils.mod

parac.mod.f90:  $(SRCDIR)/parac.mod.F90
parac.mod.o:    parac.mod.f90

para_global.mod.f90:$(SRCDIR)/para_global.mod.F90
para_global.mod.o:para_global.mod.f90

part_1d.mod.f90:$(SRCDIR)/part_1d.mod.F90
part_1d.mod.o:  part_1d.mod.f90 kinds.mod mp_interface.mod

pbc_utils.mod.f90:$(SRCDIR)/pbc_utils.mod.F90
pbc_utils.mod.o:pbc_utils.mod.f90 bc.mod clas.mod error_handling.mod \
                isos.mod kinds.mod metr.mod mm_dimmod.mod system.mod

pcgrad_driver.mod.f90:$(SRCDIR)/pcgrad_driver.mod.F90
pcgrad_driver.mod.o:pcgrad_driver.mod.f90 dotp_utils.mod elct.mod \
                ener.mod error_handling.mod forcedr_driver.mod \
                func.mod hubbardu.mod kinds.mod mm_dim_utils.mod \
                mm_dimmod.mod mm_input.mod mm_qmmm_forcedr_utils.mod \
                mp_interface.mod norm.mod nvtx_utils.mod ortho_utils.mod \
                parac.mod pslo.mod rnlsm_utils.mod rscpot_utils.mod \
                spin.mod system.mod timer.mod tpar.mod vdwcmod.mod

pcgrad_p_utils.mod.f90:$(SRCDIR)/pcgrad_p_utils.mod.F90
pcgrad_p_utils.mod.o:pcgrad_p_utils.mod.f90 csize_utils.mod \
                error_handling.mod fft_maxfft.mod fnonloc_utils.mod \
                geq0mod.mod kinds.mod mp_interface.mod norm.mod \
                parac.mod perturbation_p_utils.mod response_pmod.mod \
                rnlsm_utils.mod simple_model_p_utils.mod system.mod \
                timer.mod utils.mod vpsi_utils.mod zeroing_utils.mod \
                kinds.mod error_handling.mod fft_maxfft.mod \
                timer.mod system.mod parac.mod elct.mod response_pmod.mod \
                vpsi_utils.mod fnonloc_utils.mod rnlsm_utils.mod \
                spin.mod zeroing_utils.mod

pcgrad_utils.mod.f90:$(SRCDIR)/pcgrad_utils.mod.F90
pcgrad_utils.mod.o:pcgrad_utils.mod.f90 forcedr_utils.mod mm_dim_utils.mod \
                mm_dimmod.mod mm_input.mod ortho_utils.mod \
                rhoofr_utils.mod rnlsm_utils.mod rscpot_utils.mod \
                system.mod

pert_kpoint_p_utils.mod.f90:$(SRCDIR)/pert_kpoint_p_utils.mod.F90
pert_kpoint_p_utils.mod.o:pert_kpoint_p_utils.mod.f90 coor.mod \
                dnlpdk_p_utils.mod elct.mod ener.mod error_handling.mod \
                forces_p_utils.mod gvec.mod h0psi1_p_utils.mod \
                ions.mod kinds.mod kpert_potential_p_utils.mod \
                kpnt.mod kpts.mod ks_ener_p_utils.mod matrix_p_utils.mod \
                mp_interface.mod nlps.mod parac.mod perturbation_p_utils.mod \
                response_pmod.mod restart_p_utils.mod rhoofr_p_utils.mod \
                rhoofr_utils.mod rnl_dk_p_utils.mod rnlsm_utils.mod \
                ropt.mod rscpot_utils.mod rwfopt_p_utils.mod \
                sfac.mod spin.mod system.mod timer.mod up3_p_utils.mod \
                zeroing_utils.mod

perturbation_p_utils.mod.f90:$(SRCDIR)/perturbation_p_utils.mod.F90
perturbation_p_utils.mod.o:perturbation_p_utils.mod.f90 coor.mod \
                dotp_utils.mod error_handling.mod forces_driver.mod \
                kinds.mod mp_interface.mod ovlap_utils.mod \
                parac.mod response_pmod.mod rotate_utils.mod \
                spin.mod system.mod timer.mod zeroing_utils.mod

phfac_utils.mod.f90:$(SRCDIR)/phfac_utils.mod.F90
phfac_utils.mod.o:phfac_utils.mod.f90 cppt.mod error_handling.mod \
                gvec.mod ions.mod kinds.mod kpnt.mod kpts.mod \
                parac.mod prmem_utils.mod sfac.mod system.mod \
                timer.mod

phonons_p_utils.mod.f90:$(SRCDIR)/phonons_p_utils.mod.F90
phonons_p_utils.mod.o:phonons_p_utils.mod.f90 adat.mod coor.mod \
                cotr.mod d_mat_p_utils.mod eicalc_utils.mod \
                elct.mod error_handling.mod fnonloc_p_utils.mod \
                forces_p_utils.mod hessout_utils.mod ions.mod \
                kinds.mod kpnt.mod mp_interface.mod nlps.mod \
                parac.mod response_pmod.mod rhoofr_p_utils.mod \
                rhoofr_utils.mod rmas.mod rnlsm_p_utils.mod \
                rnlsm_utils.mod ropt.mod rscpot_utils.mod rwfopt_p_utils.mod \
                secder_utils.mod sfac.mod spin.mod symm.mod \
                symtrz_utils.mod system.mod timer.mod utils.mod \
                zeroing_utils.mod

pi_cntl_utils.mod.f90:$(SRCDIR)/pi_cntl_utils.mod.F90
pi_cntl_utils.mod.o:pi_cntl_utils.mod.f90 cotr.mod error_handling.mod \
                inscan_utils.mod ions.mod kinds.mod parac.mod \
                pimd.mod readsr_utils.mod store_types.mod system.mod

pi_diag_utils.mod.f90:$(SRCDIR)/pi_diag_utils.mod.F90
pi_diag_utils.mod.o:pi_diag_utils.mod.f90 andp.mod anneal_utils.mod \
                atwf.mod calc_alm_utils.mod cnst.mod cnstpr_utils.mod \
                comvel_utils.mod comvelmod.mod coor.mod copot_utils.mod \
                cotr.mod detdof_utils.mod dispp_utils.mod dynit_utils.mod \
                ekinpp_utils.mod elct.mod ener.mod error_handling.mod \
                evirial_utils.mod fharm_utils.mod fileopen_utils.mod \
                fileopenmod.mod filnmod.mod fint.mod forcep_utils.mod \
                forces_diag_utils.mod freqs_utils.mod geofile_utils.mod \
                getcor_utils.mod getfnm_utils.mod getfu_utils.mod \
                getgyr_utils.mod global_utils.mod ions.mod \
                kinds.mod kpts.mod localize_utils.mod machine.mod \
                md_driver.mod mm_extrap.mod moverho_utils.mod \
                mp_interface.mod nlcc.mod nose.mod noseng_utils.mod \
                nosepa_utils.mod noseup_utils.mod nospinit_utils.mod \
                parac.mod phfac_utils.mod pi_md_utils.mod pimd.mod \
                pinmtrans_utils.mod pi_stress_utils.mod poin.mod \
                posupi_utils.mod printave_utils.mod printp_utils.mod \
                prmem_utils.mod prtgyr_utils.mod pslo.mod puttau_utils.mod \
                ranp_utils.mod rattle_utils.mod readsr_utils.mod \
                reshaper.mod rhoofr_c_utils.mod rhoofr_utils.mod \
                rinitwf_driver.mod rinvel_utils.mod rnlsm_utils.mod \
                ropt.mod rotvel_utils.mod rscvp_utils.mod rv30_utils.mod \
                setirec_utils.mod shake_utils.mod soft.mod \
                spin.mod stagetrans_utils.mod store_types.mod \
                strs.mod system.mod testex_utils.mod teststore_utils.mod \
                totstr_utils.mod vdw_utils.mod vdwcmod.mod \
                velupi_utils.mod wr_temps_utils.mod wrener_utils.mod \
                wrgeo_utils.mod wv30_utils.mod zeroing_utils.mod

pi_init_utils.mod.f90:$(SRCDIR)/pi_init_utils.mod.F90
pi_init_utils.mod.o:pi_init_utils.mod.f90 adat.mod cnst.mod \
                error_handling.mod ions.mod kinds.mod mp_multiple_comm_init.mod \
                parac.mod pimd.mod rmas.mod system.mod timer.mod \
                utils.mod wann.mod zeroing_utils.mod

pimd.mod.f90:   $(SRCDIR)/pimd.mod.F90
pimd.mod.o:     pimd.mod.f90 kinds.mod system.mod

pi_mdpt_utils.mod.f90:$(SRCDIR)/pi_mdpt_utils.mod.F90
pi_mdpt_utils.mod.o:pi_mdpt_utils.mod.f90 atwf.mod ddip.mod \
                elct.mod error_handling.mod fnlalloc_utils.mod \
                kinds.mod kpts.mod mp_interface.mod parac.mod \
                pi_diag_utils.mod pi_md_utils.mod pimd.mod \
                prmem_utils.mod pslo.mod system.mod testex_utils.mod \
                timer.mod utils.mod vdwcmod.mod zeroing_utils.mod

pi_md_utils.mod.f90:$(SRCDIR)/pi_md_utils.mod.F90
pi_md_utils.mod.o:pi_md_utils.mod.f90 anneal_utils.mod cnst.mod \
                cnstpr_utils.mod comvel_utils.mod comvelmod.mod \
                coor.mod copot_utils.mod cotr.mod ddipo_utils.mod \
                deort_utils.mod detdof_utils.mod dispp_utils.mod \
                dynit_utils.mod ekinpp_utils.mod elct.mod ener.mod \
                error_handling.mod evirial_utils.mod fharm_utils.mod \
                fileopen_utils.mod fileopenmod.mod filnmod.mod \
                forcedr_driver.mod forcep_utils.mod freqs_utils.mod \
                geofile_utils.mod getcor_utils.mod getfnm_utils.mod \
                getfu_utils.mod getgyr_utils.mod global_utils.mod \
                ions.mod kinds.mod localize_utils.mod machine.mod \
                mdmain_utils.mod mp_interface.mod nlcc.mod \
                nose.mod noseinit_utils.mod noseng_utils.mod \
                nosepa_utils.mod noseup_utils.mod nospinit_utils.mod \
                parac.mod phfac_utils.mod pimd.mod pinmtrans_utils.mod \
                pi_stress_utils.mod posupa_utils.mod posupi_utils.mod \
                printave_utils.mod printp_utils.mod prmem_utils.mod \
                prtgyr_utils.mod pslo.mod puttau_utils.mod \
                ranp_utils.mod rattle_utils.mod readsr_utils.mod \
                rekine_utils.mod reshaper.mod rhopri_utils.mod \
                rinvel_utils.mod ropt.mod rortv_utils.mod rotvel_utils.mod \
                rscve_utils.mod rscvp_utils.mod rv30_utils.mod \
                setirec_utils.mod shake_utils.mod soft.mod \
                spin.mod stagetrans_utils.mod store_types.mod \
                strs.mod system.mod testex_utils.mod teststore_utils.mod \
                timer.mod totstr_utils.mod vdw_utils.mod vdwcmod.mod \
                velupa_utils.mod velupi_utils.mod wannier_print_utils.mod \
                wr_temps_utils.mod wrener_utils.mod wrgeo_utils.mod \
                wv30_utils.mod zeroing_utils.mod

pimd_utils.mod.f90:$(SRCDIR)/pimd_utils.mod.F90
pimd_utils.mod.o:pimd_utils.mod.f90 error_handling.mod pimd.mod

pinmtrans_utils.mod.f90:$(SRCDIR)/pinmtrans_utils.mod.F90
pinmtrans_utils.mod.o:pinmtrans_utils.mod.f90 ions.mod kinds.mod \
                pimd.mod zeroing_utils.mod

pi_npt_bomd_utils.mod.f90:$(SRCDIR)/pi_npt_bomd_utils.mod.F90
pi_npt_bomd_utils.mod.o:pi_npt_bomd_utils.mod.f90 andp.mod \
                anneal_utils.mod atwf.mod calc_alm_utils.mod \
                cnst.mod cnstpr_utils.mod comvel_utils.mod \
                comvelmod.mod coor.mod copot_utils.mod cotr.mod \
                detdof_utils.mod dispp_utils.mod dynit_utils.mod \
                ekinpp_utils.mod elct.mod ener.mod error_handling.mod \
                evirial_utils.mod fharm_utils.mod fileopen_utils.mod \
                fileopenmod.mod filnmod.mod fint.mod forcep_utils.mod \
                forces_diag_utils.mod freqs_utils.mod geofile_utils.mod \
                getcor_utils.mod getfnm_utils.mod getfu_utils.mod \
                getgyr_utils.mod global_utils.mod ions.mod \
                kinds.mod kpts.mod localize_utils.mod machine.mod \
                metr.mod mm_extrap.mod moverho_utils.mod mp_interface.mod \
                newcell_utils.mod nlcc.mod noscinit_utils.mod \
                nose.mod noseng_utils.mod nosepa_utils.mod \
                noseup_utils.mod nospinit_utils.mod parac.mod \
                phfac_utils.mod pi_md_utils.mod pimd.mod pinmtrans_utils.mod \
                pi_stress_utils.mod poin.mod posupi_utils.mod \
                prbomd_utils.mod prcp.mod printave_utils.mod \
                printp_utils.mod prmem_utils.mod prtgyr_utils.mod \
                pslo.mod puttau_utils.mod ranp_utils.mod rattle_utils.mod \
                readsr_utils.mod reshaper.mod rhoofr_c_utils.mod \
                rhoofr_utils.mod rinitwf_driver.mod rinvel_utils.mod \
                rnlsm_utils.mod ropt.mod rotvel_utils.mod rscvp_utils.mod \
                rv30_utils.mod setirec_utils.mod setsc_utils.mod \
                shake_utils.mod shock.mod soft.mod spin.mod \
                stagetrans_utils.mod store_types.mod strs.mod \
                system.mod testex_utils.mod teststore_utils.mod \
                totstr_utils.mod vdw_utils.mod vdwcmod.mod \
                velupi_utils.mod wr_temps_utils.mod wrener_utils.mod \
                wrgeo_utils.mod wv30_utils.mod zeroing_utils.mod

pi_npt_cpmd_utils.mod.f90:$(SRCDIR)/pi_npt_cpmd_utils.mod.F90
pi_npt_cpmd_utils.mod.o:pi_npt_cpmd_utils.mod.f90 anneal_utils.mod \
                cnst.mod cnstpr_utils.mod comvel_utils.mod \
                comvelmod.mod coor.mod copot_utils.mod cotr.mod \
                ddipo_utils.mod deort_utils.mod detdof_utils.mod \
                dispp_utils.mod dynit_utils.mod ekinpp_utils.mod \
                elct.mod ener.mod error_handling.mod evirial_utils.mod \
                fharm_utils.mod fileopen_utils.mod fileopenmod.mod \
                filnmod.mod forcedr_driver.mod forcep_utils.mod \
                freqs_utils.mod geofile_utils.mod getcor_utils.mod \
                getfnm_utils.mod getfu_utils.mod getgyr_utils.mod \
                global_utils.mod ions.mod kinds.mod localize_utils.mod \
                machine.mod metr.mod mp_interface.mod newcell_utils.mod \
                nlcc.mod nose.mod noscinit_utils.mod noseinit_utils.mod \
                noseng_utils.mod nosepa_utils.mod noseup_utils.mod \
                nospinit_utils.mod parac.mod phfac_utils.mod \
                pimd.mod pi_md_utils.mod pinmtrans_utils.mod \
                pi_stress_utils.mod posupa_utils.mod posupi_utils.mod \
                prcp.mod prcpmd_utils.mod printave_utils.mod \
                printp_utils.mod prmem_utils.mod prtgyr_utils.mod \
                pslo.mod puttau_utils.mod ranp_utils.mod rattle_utils.mod \
                readsr_utils.mod rekine_utils.mod reshaper.mod \
                rhopri_utils.mod rinvel_utils.mod ropt.mod \
                rortv_utils.mod rotvel_utils.mod rscve_utils.mod \
                rscvp_utils.mod rv30_utils.mod setirec_utils.mod \
                setsc_utils.mod shake_utils.mod shock.mod soft.mod \
                spin.mod stagetrans_utils.mod store_types.mod \
                strs.mod system.mod testex_utils.mod teststore_utils.mod \
                timer.mod totstr_utils.mod vdw_utils.mod vdwcmod.mod \
                velupa_utils.mod velupi_utils.mod wannier_print_utils.mod \
                wr_temps_utils.mod wrener_utils.mod wrgeo_utils.mod \
                wv30_utils.mod zeroing_utils.mod

pi_prpt_utils.mod.f90:$(SRCDIR)/pi_prpt_utils.mod.F90
pi_prpt_utils.mod.o:pi_prpt_utils.mod.f90 atwf.mod ddip.mod \
                elct.mod error_handling.mod fnlalloc_utils.mod \
                kinds.mod kpts.mod mp_interface.mod parac.mod \
                pimd.mod pi_npt_bomd_utils.mod pi_npt_cpmd_utils.mod \
                prmem_utils.mod pslo.mod system.mod testex_utils.mod \
                timer.mod utils.mod vdwcmod.mod zeroing_utils.mod

pi_stress_utils.mod.f90:$(SRCDIR)/pi_stress_utils.mod.F90
pi_stress_utils.mod.o:pi_stress_utils.mod.f90 cnst.mod fileopen_utils.mod \
                fileopenmod.mod ions.mod kinds.mod metr.mod \
                mp_interface.mod parac.mod pbc_utils.mod pimd.mod \
                prcp.mod rmas.mod ropt.mod store_types.mod \
                symtrz_utils.mod system.mod zeroing_utils.mod

pi_wf_utils.mod.f90:$(SRCDIR)/pi_wf_utils.mod.F90
pi_wf_utils.mod.o:pi_wf_utils.mod.f90 coor.mod error_handling.mod \
                filnmod.mod parac.mod pimd.mod readsr_utils.mod \
                repgen_utils.mod rreadf_utils.mod store_types.mod \
                system.mod timer.mod wfopts_utils.mod

plotband.f90:   $(SRCDIR)/plotband.F90
plotband.o:     plotband.f90 inscan_utils.mod machine.mod

pm_cntl_utils.mod.f90:$(SRCDIR)/pm_cntl_utils.mod.F90
pm_cntl_utils.mod.o:pm_cntl_utils.mod.f90 cotr.mod error_handling.mod \
                inscan_utils.mod ions.mod kinds.mod mfep.mod \
                parac.mod pimd.mod readsr_utils.mod store_types.mod \
                system.mod

pm_gmopts_utils.mod.f90:$(SRCDIR)/pm_gmopts_utils.mod.F90
pm_gmopts_utils.mod.o:pm_gmopts_utils.mod.f90 atwf.mod coor.mod \
                elct.mod error_handling.mod filnmod.mod fnlalloc_utils.mod \
                kinds.mod linres.mod lr_in_utils.mod parac.mod \
                pimd.mod prmem_utils.mod pslo.mod readsr_utils.mod \
                rreadf_utils.mod store_types.mod system.mod \
                timer.mod zeroing_utils.mod

pm_init_utils.mod.f90:$(SRCDIR)/pm_init_utils.mod.F90
pm_init_utils.mod.o:pm_init_utils.mod.f90 cotr.mod envj.mod \
                error_handling.mod mfep.mod mp_interface.mod \
                mp_multiple_comm_init.mod parac.mod pimd.mod \
                timer.mod wann.mod

pm_mdpt_utils.mod.f90:$(SRCDIR)/pm_mdpt_utils.mod.F90
pm_mdpt_utils.mod.o:pm_mdpt_utils.mod.f90 atwf.mod bsym.mod \
                cl_init_utils.mod clas.mod cnst.mod coninp_utils.mod \
                coor.mod cotr.mod ddip.mod elct.mod error_handling.mod \
                fileopen_utils.mod fileopenmod.mod filnmod.mod \
                fnlalloc_utils.mod fusion_utils.mod jrotation_utils.mod \
                kinds.mod kpts.mod linres.mod lr_in_utils.mod \
                md_driver.mod mdclas_utils.mod mdfile_utils.mod \
                mdmain_utils.mod mdshop_bo_utils.mod mdshop_cp_utils.mod \
                meta_multiple_walkers_utils.mod mfep.mod mm_dim_utils.mod \
                mm_dimmod.mod mm_mddiag_utils.mod mm_mdmain_utils.mod \
                mm_mdshop_bo_utils.mod mm_mdshop_cp_utils.mod \
                mp_interface.mod parac.mod pimd.mod printpmod.mod \
                prmem_utils.mod pslo.mod readsr_utils.mod rreadf_utils.mod \
                soft.mod store_types.mod system.mod timer.mod \
                utils.mod vdwcmod.mod zeroing_utils.mod

pm_wf_utils.mod.f90:$(SRCDIR)/pm_wf_utils.mod.F90
pm_wf_utils.mod.o:pm_wf_utils.mod.f90 coor.mod error_handling.mod \
                filnmod.mod parac.mod pimd.mod readsr_utils.mod \
                rreadf_utils.mod store_types.mod system.mod \
                timer.mod wfopts_utils.mod

pnosmove_utils.mod.f90:$(SRCDIR)/pnosmove_utils.mod.F90
pnosmove_utils.mod.o:pnosmove_utils.mod.f90 cnst.mod ions.mod \
                kinds.mod nose.mod system.mod

poin.mod.f90:   $(SRCDIR)/poin.mod.F90
poin.mod.o:     poin.mod.f90 kinds.mod

pola.mod.f90:   $(SRCDIR)/pola.mod.F90
pola.mod.o:     pola.mod.f90 kinds.mod

polarise_utils.mod.f90:$(SRCDIR)/polarise_utils.mod.F90
polarise_utils.mod.o:polarise_utils.mod.f90 atimesmod.mod calc_pij_utils.mod \
                clinbcg_utils.mod cnst.mod cppt.mod dotp_utils.mod \
                error_handling.mod fft_maxfft.mod fftmain_utils.mod \
                fftnew_utils.mod geq0mod.mod ions.mod kinds.mod \
                kpnt.mod kpts.mod machine.mod parac.mod pola.mod \
                projv_utils.mod sfac.mod spin.mod system.mod \
                timer.mod utils.mod wrener_utils.mod zeroing_utils.mod

posupa_utils.mod.f90:$(SRCDIR)/posupa_utils.mod.F90
posupa_utils.mod.o:posupa_utils.mod.f90 dotp_utils.mod error_handling.mod \
                harm.mod jrotation_utils.mod kinds.mod kpts.mod \
                mp_interface.mod parac.mod pslo.mod ropt.mod \
                rortog_utils.mod rotate_utils.mod spin.mod \
                system.mod timer.mod tpar.mod

posupi_utils.mod.f90:$(SRCDIR)/posupi_utils.mod.F90
posupi_utils.mod.o:posupi_utils.mod.f90 cnst.mod ions.mod jacobi_utils.mod \
                kinds.mod metr.mod mp_interface.mod parac.mod \
                readsr_utils.mod store_types.mod system.mod \
                tpar.mod zeroing_utils.mod

potfor_utils.mod.f90:$(SRCDIR)/potfor_utils.mod.F90
potfor_utils.mod.o:potfor_utils.mod.f90 cppt.mod geq0mod.mod \
                ions.mod kinds.mod sfac.mod system.mod timer.mod

potmed_utils.mod.f90:$(SRCDIR)/potmed_utils.mod.F90
potmed_utils.mod.o:potmed_utils.mod.f90 adat.mod atomc_utils.mod \
                ions.mod kinds.mod metr.mod mp_interface.mod \
                parac.mod pbc_utils.mod system.mod timer.mod

ppener_utils.mod.f90:$(SRCDIR)/ppener_utils.mod.F90
ppener_utils.mod.o:ppener_utils.mod.f90 cppt.mod geq0mod.mod \
                kinds.mod nvtx_utils.mod simulmod.mod system.mod \
                timer.mod

prbomd_utils.mod.f90:$(SRCDIR)/prbomd_utils.mod.F90
prbomd_utils.mod.o:prbomd_utils.mod.f90 andp.mod andr.mod anneal_utils.mod \
                atwf.mod bsym.mod calc_alm_utils.mod cnst.mod \
                cnst_dyn.mod comvel_utils.mod comvelmod.mod \
                coor.mod copot_utils.mod cotr.mod csize_utils.mod \
                ddipo_utils.mod detdof_utils.mod dispp_utils.mod \
                dynit_utils.mod ekinpp_utils.mod elct.mod ener.mod \
                error_handling.mod extrap_utils.mod fileopen_utils.mod \
                fileopenmod.mod finalp_utils.mod fint.mod forcep_utils.mod \
                forces_diag_utils.mod geofile_utils.mod gsize_utils.mod \
                initrun_driver.mod initrun_utils.mod ions.mod \
                kinds.mod kpts.mod localize_utils.mod machine.mod \
                meta_cell_utils.mod meta_colvar_inp_utils.mod \
                meta_exlagr_methods.mod meta_exlagr_utils.mod \
                metr.mod mm_extrap.mod moverho_utils.mod mp_interface.mod \
                newcell_utils.mod nlcc.mod norm.mod nose.mod \
                parac.mod phfac_utils.mod poin.mod posupi_utils.mod \
                prcp.mod printave_utils.mod printp_utils.mod \
                proja_utils.mod puttau_utils.mod rattle_utils.mod \
                resetac_utils.mod rhopri_utils.mod rinvel_utils.mod \
                ropt.mod rscvp_utils.mod setirec_utils.mod \
                setsc_utils.mod shake_utils.mod soft.mod spin.mod \
                store_types.mod str2.mod system.mod testex_utils.mod \
                teststore_utils.mod timer.mod totstr_utils.mod \
                tpar.mod vdwcmod.mod velupi_utils.mod wannier_print_utils.mod \
                wrener_utils.mod wrgeo_utils.mod wv30_utils.mod \
                zeroing_utils.mod

prcnosmove_utils.mod.f90:$(SRCDIR)/prcnosmove_utils.mod.F90
prcnosmove_utils.mod.o:prcnosmove_utils.mod.f90 cnst.mod ions.mod \
                jacobi_utils.mod kinds.mod metr.mod nose.mod \
                system.mod zeroing_utils.mod

prcpmd_utils.mod.f90:$(SRCDIR)/prcpmd_utils.mod.F90
prcpmd_utils.mod.o:prcpmd_utils.mod.f90 anneal_utils.mod cnst.mod \
                cnst_dyn.mod comvel_utils.mod comvelmod.mod \
                coor.mod copot_utils.mod cotr.mod csize_utils.mod \
                ddipo_utils.mod deort_utils.mod detdof_utils.mod \
                dispp_utils.mod dynit_utils.mod ekinpp_utils.mod \
                elct.mod ener.mod error_handling.mod fileopen_utils.mod \
                fileopenmod.mod finalp_utils.mod forcedr_driver.mod \
                forcedr_utils.mod forcep_utils.mod freqs_utils.mod \
                geofile_utils.mod geq0mod.mod gsize_utils.mod \
                initrun_driver.mod initrun_utils.mod ions.mod \
                kinds.mod kpts.mod localize_utils.mod machine.mod \
                meta_cell_utils.mod meta_colvar_inp_utils.mod \
                meta_colvar_utils.mod meta_exlagr_methods.mod \
                meta_exlagr_utils.mod metr.mod mp_interface.mod \
                newcell_utils.mod nlcc.mod norm.mod nose.mod \
                ortho_utils.mod parac.mod phfac_utils.mod posupa_utils.mod \
                posupi_utils.mod prcp.mod printave_utils.mod \
                printp_utils.mod proja_utils.mod pslo.mod puttau_utils.mod \
                quenbo_utils.mod rattle_utils.mod rekine_utils.mod \
                resetac_utils.mod rhopri_utils.mod rinvel_utils.mod \
                ropt.mod rortv_utils.mod rscve_utils.mod rscvp_utils.mod \
                setirec_utils.mod setsc_utils.mod shake_utils.mod \
                soft.mod spin.mod store_types.mod str2.mod \
                system.mod testex_utils.mod teststore_utils.mod \
                totstr_utils.mod tpar.mod utils.mod vdwcmod.mod \
                velupa_utils.mod velupi_utils.mod wannier_print_utils.mod \
                wrener_utils.mod wrgeo_utils.mod wv30_utils.mod \
                zeroing_utils.mod

prcp.mod.f90:   $(SRCDIR)/prcp.mod.F90
prcp.mod.o:     prcp.mod.f90 kinds.mod

prden.mod.f90:  $(SRCDIR)/prden.mod.F90
prden.mod.o:    prden.mod.f90 kinds.mod

prep_forcematch_utils.mod.f90:$(SRCDIR)/prep_forcematch_utils.mod.F90
prep_forcematch_utils.mod.o:prep_forcematch_utils.mod.f90 bsym.mod \
                ddip.mod elct.mod error_handling.mod fnlalloc_utils.mod \
                jrotation_utils.mod kinds.mod linres.mod mm_dim_utils.mod \
                mm_dimmod.mod mm_forcematch_utils.mod parac.mod \
                prmem_utils.mod pslo.mod shop_rest_2.mod system.mod \
                utils.mod zeroing_utils.mod

printave_utils.mod.f90:$(SRCDIR)/printave_utils.mod.F90
printave_utils.mod.o:printave_utils.mod.f90 kinds.mod parac.mod \
                system.mod

printfor_utils.mod.f90:$(SRCDIR)/printfor_utils.mod.F90
printfor_utils.mod.o:printfor_utils.mod.f90 fileopen_utils.mod \
                fileopenmod.mod ions.mod kinds.mod parac.mod \
                ropt.mod store_types.mod system.mod

printp.mod.f90: $(SRCDIR)/printp.mod.F90
printp.mod.o:   printp.mod.f90

printp_utils.mod.f90:$(SRCDIR)/printp_utils.mod.F90
printp_utils.mod.o:printp_utils.mod.f90 adat.mod cell.mod clas.mod \
                cnst.mod coninp_utils.mod cotr.mod ddip.mod \
                error_handling.mod fileopen_utils.mod fileopenmod.mod \
                ions.mod kinds.mod meta_multiple_walkers_utils.mod \
                metr.mod mfep.mod mm_dim_utils.mod mm_dimmod.mod \
                movi.mod mw.mod parac.mod pimd.mod printpmod.mod \
                readsr_utils.mod rmas.mod ropt.mod store_types.mod \
                strs.mod system.mod bicanonicalCpmd.mod

prmdfile_utils.mod.f90:$(SRCDIR)/prmdfile_utils.mod.F90
prmdfile_utils.mod.o:prmdfile_utils.mod.f90 andp.mod andr.mod \
                anneal_utils.mod atwf.mod calc_alm_utils.mod \
                cnst.mod comvel_utils.mod comvelmod.mod coor.mod \
                copot_utils.mod cotr.mod ddipo_utils.mod detdof_utils.mod \
                dispp_utils.mod do_perturbation_p_utils.mod \
                dynit_utils.mod ekinpp_utils.mod elct.mod ener.mod \
                error_handling.mod extrap_utils.mod fileopen_utils.mod \
                fileopenmod.mod finalp_utils.mod fint.mod forcep_utils.mod \
                forces_diag_utils.mod geofile_utils.mod gsize_utils.mod \
                initrun_driver.mod initrun_utils.mod kinds.mod \
                kpts.mod linres.mod localize_utils.mod lr_tddft_utils.mod \
                machine.mod metr.mod mm_extrap.mod moverho_utils.mod \
                mp_interface.mod newcell_utils.mod nlcc.mod \
                norm.mod nose.mod parac.mod phfac_utils.mod \
                poin.mod posupi_utils.mod prcp.mod printave_utils.mod \
                printfor_utils.mod printp_utils.mod proppt_utils.mod \
                puttau_utils.mod rattle_utils.mod resetac_utils.mod \
                rhopri_utils.mod rinvel_utils.mod ropt.mod \
                rscvp_utils.mod sample_utils.mod setirec_utils.mod \
                setsc_utils.mod shake_utils.mod soft.mod spin.mod \
                store_types.mod system.mod testex_utils.mod \
                teststore_utils.mod totstr_utils.mod vdwcmod.mod \
                velupi_utils.mod wrener_utils.mod wv30_utils.mod \
                zeroing_utils.mod

prmem_utils.mod.f90:$(SRCDIR)/prmem_utils.mod.F90
prmem_utils.mod.o:prmem_utils.mod.f90 envj.mod kinds.mod machine.mod \
                parac.mod

prng.mod.f90:   $(SRCDIR)/prng.mod.F90
prng.mod.o:     prng.mod.f90 kinds.mod

prng_utils.mod.f90:$(SRCDIR)/prng_utils.mod.F90
prng_utils.mod.o:prng_utils.mod.f90 kinds.mod parac.mod prng.mod \
                system.mod

proja_utils.mod.f90:$(SRCDIR)/proja_utils.mod.F90
proja_utils.mod.o:proja_utils.mod.f90 csize_utils.mod dotp_utils.mod \
                error_handling.mod kinds.mod mp_interface.mod \
                norm.mod ovlap_utils.mod parac.mod pslo.mod \
                rotate_utils.mod spin.mod spsi_utils.mod system.mod \
                timer.mod

projv_utils.mod.f90:$(SRCDIR)/projv_utils.mod.F90
projv_utils.mod.o:projv_utils.mod.f90 dotp_utils.mod kinds.mod \
                kpts.mod pola.mod system.mod

propin_utils.mod.f90:$(SRCDIR)/propin_utils.mod.F90
propin_utils.mod.o:propin_utils.mod.f90 cnst.mod condu.mod \
                cores.mod error_handling.mod g_loc.mod inscan_utils.mod \
                ldosmod.mod lodp.mod mp_interface.mod parac.mod \
                pola.mod prop.mod readsr_utils.mod system.mod \
                wann.mod zeroing_utils.mod

prop.mod.f90:   $(SRCDIR)/prop.mod.F90
prop.mod.o:     prop.mod.f90 kinds.mod system.mod

proppt_utils.mod.f90:$(SRCDIR)/proppt_utils.mod.F90
proppt_utils.mod.o:proppt_utils.mod.f90 adat.mod atomc_utils.mod \
                cnst.mod condu.mod conduct_utils.mod coor.mod \
                core_spect_utils.mod cores.mod cppt.mod ddip.mod \
                ddipo_utils.mod difrho_utils.mod dipo_utils.mod \
                dipomod.mod dist_prowfn_utils.mod eicalc_utils.mod \
                elct.mod elstpo_utils.mod ener.mod error_handling.mod \
                espchg_utils.mod exdipo_utils.mod fft_maxfft.mod \
                fftmain_utils.mod fftnew_utils.mod fileopen_utils.mod \
                fileopenmod.mod fnlalloc_utils.mod forcedr_driver.mod \
                forcep_utils.mod forces_utils.mod g_loc.mod \
                g_loc_dr_utils.mod geq0mod.mod hip_utils.mod \
                ions.mod isos.mod kddipo_utils.mod kinds.mod \
                kpts.mod ldos_utils.mod ldosmod.mod localize_utils.mod \
                lodipo_utils.mod lodp.mod mm_dim_utils.mod \
                mm_dimmod.mod molorb_utils.mod mp_interface.mod \
                ortho_utils.mod parac.mod perturbation_p_utils.mod \
                phfac_utils.mod poin.mod pola.mod polarise_utils.mod \
                potmed_utils.mod prmem_utils.mod prop.mod propin_utils.mod \
                prowfn_utils.mod proylm_utils.mod pslo.mod \
                readsr_utils.mod rho1ofr_utils.mod rhoofr_c_utils.mod \
                rhoofr_utils.mod rnlsm_utils.mod ropt.mod rscpot_utils.mod \
                rv30_utils.mod setirec_utils.mod soft.mod spin.mod \
                store_types.mod system.mod testex_utils.mod \
                utils.mod vdwcmod.mod vofrho_utils.mod wann.mod \
                wannier_print_utils.mod wc_dos_utils.mod zeroing_utils.mod

prowfn_utils.mod.f90:$(SRCDIR)/prowfn_utils.mod.F90
prowfn_utils.mod.o:prowfn_utils.mod.f90 adat.mod atwf.mod augchg_utils.mod \
                cmaos_utils.mod dotp_utils.mod elct.mod error_handling.mod \
                fileopen_utils.mod fileopenmod.mod forcep_utils.mod \
                ions.mod kinds.mod machine.mod mp_interface.mod \
                ovlap_utils.mod parac.mod prden.mod prmem_utils.mod \
                prop.mod pslo.mod rnlsm_utils.mod setbasis_utils.mod \
                sfac.mod spin.mod summat_utils.mod system.mod \
                timer.mod utils.mod zeroing_utils.mod

proylm_utils.mod.f90:$(SRCDIR)/proylm_utils.mod.F90
proylm_utils.mod.o:proylm_utils.mod.f90 bessm_utils.mod cnst.mod \
                cppt.mod dotp_utils.mod error_handling.mod \
                fileopen_utils.mod fileopenmod.mod gvec.mod \
                kinds.mod mp_interface.mod parac.mod prop.mod \
                system.mod timer.mod zeroing_utils.mod

prpcmove_utils.mod.f90:$(SRCDIR)/prpcmove_utils.mod.F90
prpcmove_utils.mod.o:prpcmove_utils.mod.f90 ions.mod jacobi_utils.mod \
                kinds.mod metr.mod system.mod zeroing_utils.mod

prpcnosmove_utils.mod.f90:$(SRCDIR)/prpcnosmove_utils.mod.F90
prpcnosmove_utils.mod.o:prpcnosmove_utils.mod.f90 cnst.mod \
                ions.mod jacobi_utils.mod kinds.mod metr.mod \
                nose.mod system.mod zeroing_utils.mod

prpnosmove_utils.mod.f90:$(SRCDIR)/prpnosmove_utils.mod.F90
prpnosmove_utils.mod.o:prpnosmove_utils.mod.f90 cnst.mod ions.mod \
                jacobi_utils.mod kinds.mod metr.mod nose.mod \
                system.mod zeroing_utils.mod

prpt_utils.mod.f90:$(SRCDIR)/prpt_utils.mod.F90
prpt_utils.mod.o:prpt_utils.mod.f90 atwf.mod ddip.mod elct.mod \
                error_handling.mod fnlalloc_utils.mod jrotation_utils.mod \
                kinds.mod npt_md_utils.mod parac.mod prbomd_utils.mod \
                prcpmd_utils.mod prmdfile_utils.mod prmem_utils.mod \
                pslo.mod system.mod utils.mod vdwcmod.mod zeroing_utils.mod

prtgyr_utils.mod.f90:$(SRCDIR)/prtgyr_utils.mod.F90
prtgyr_utils.mod.o:prtgyr_utils.mod.f90 adat.mod cnst.mod ions.mod \
                parac.mod pimd.mod

pslo.mod.f90:   $(SRCDIR)/pslo.mod.F90
pslo.mod.o:     pslo.mod.f90 system.mod

pstat.mod.f90:  $(SRCDIR)/pstat.mod.F90
pstat.mod.o:    pstat.mod.f90 kinds.mod

ptheory_utils.mod.f90:$(SRCDIR)/ptheory_utils.mod.F90
ptheory_utils.mod.o:ptheory_utils.mod.f90 error_handling.mod \
                fint.mod hpsi_utils.mod kinds.mod kpts.mod \
                rgs_utils.mod system.mod timer.mod

purge_utils.mod.f90:$(SRCDIR)/purge_utils.mod.F90
purge_utils.mod.o:purge_utils.mod.f90 cotr.mod error_handling.mod \
                ions.mod kinds.mod mp_interface.mod parac.mod \
                utils.mod zeroing_utils.mod

putbet_utils.mod.f90:$(SRCDIR)/putbet_utils.mod.F90
putbet_utils.mod.o:putbet_utils.mod.f90 cppt.mod dpot.mod dylmr_utils.mod \
                fitpack_utils.mod geq0mod.mod ions.mod kinds.mod \
                kpnt.mod kpts.mod nlps.mod pslo.mod qspl.mod \
                sgpp.mod str2.mod system.mod vdbp.mod ylmr_utils.mod

puttau_utils.mod.f90:$(SRCDIR)/puttau_utils.mod.F90
puttau_utils.mod.o:puttau_utils.mod.f90 cotr.mod ions.mod kinds.mod

pw_hfx_input_cnst.mod.f90:$(SRCDIR)/pw_hfx_input_cnst.mod.F90
pw_hfx_input_cnst.mod.o:pw_hfx_input_cnst.mod.f90

pw_hfx.mod.f90: $(SRCDIR)/pw_hfx.mod.F90
pw_hfx.mod.o:   pw_hfx.mod.f90 cnst.mod cp_grp_utils.mod cppt.mod \
                dotp_utils.mod error_handling.mod fft.mod fft_maxfft.mod \
                fftmain_utils.mod fftnew_utils.mod func.mod \
                geq0mod.mod hfx_utils.mod hfxmod.mod kinds.mod \
                kpts.mod machine.mod min_heap.mod mp_interface.mod \
                parac.mod part_1d.mod pw_hfx_input_cnst.mod \
                ropt.mod spin.mod state_utils.mod system.mod \
                timer.mod zeroing_utils.mod

pw_hfx_resp.mod.f90:$(SRCDIR)/pw_hfx_resp.mod.F90
pw_hfx_resp.mod.o:pw_hfx_resp.mod.f90 cnst.mod cp_grp_utils.mod \
                cppt.mod elct.mod error_handling.mod fft.mod \
                fft_maxfft.mod fftmain_utils.mod fftnew_utils.mod \
                func.mod geq0mod.mod hfxmod.mod kinds.mod kpts.mod \
                mp_interface.mod parac.mod pslo.mod pw_hfx.mod \
                pw_hfx_input_cnst.mod pw_hfx_resp_types.mod \
                pw_hfx_resp_utils.mod spin.mod state_utils.mod \
                system.mod timer.mod zeroing_utils.mod

pw_hfx_resp_types.mod.f90:$(SRCDIR)/pw_hfx_resp_types.mod.F90
pw_hfx_resp_types.mod.o:pw_hfx_resp_types.mod.f90 kinds.mod

pw_hfx_resp_utils.mod.f90:$(SRCDIR)/pw_hfx_resp_utils.mod.F90
pw_hfx_resp_utils.mod.o:pw_hfx_resp_utils.mod.f90 error_handling.mod \
                pw_hfx_resp_types.mod

qrada_s_utils.mod.f90:$(SRCDIR)/qrada_s_utils.mod.F90
qrada_s_utils.mod.o:qrada_s_utils.mod.f90 cnst.mod cppt.mod \
                dylmr_utils.mod error_handling.mod fitpack_utils.mod \
                ions.mod kinds.mod mp_interface.mod parac.mod \
                pslo.mod qspl.mod radin_utils.mod str2.mod \
                system.mod timer.mod vdbp.mod zeroing_utils.mod

qspl.mod.f90:   $(SRCDIR)/qspl.mod.F90
qspl.mod.o:     qspl.mod.f90 kinds.mod

quenbo_utils.mod.f90:$(SRCDIR)/quenbo_utils.mod.F90
quenbo_utils.mod.o:quenbo_utils.mod.f90 coor.mod elct.mod ener.mod \
                error_handling.mod kinds.mod machine.mod mm_input.mod \
                mp_interface.mod norm.mod parac.mod prmem_utils.mod \
                ropt.mod system.mod timer.mod updwf_utils.mod \
                wrener_utils.mod

qvan1_utils.mod.f90:$(SRCDIR)/qvan1_utils.mod.F90
qvan1_utils.mod.o:qvan1_utils.mod.f90 aavan.mod cnst.mod cvan.mod \
                kinds.mod nlps.mod system.mod

qvan2_utils.mod.f90:$(SRCDIR)/qvan2_utils.mod.F90
qvan2_utils.mod.o:qvan2_utils.mod.f90 aavan.mod cppt.mod cvan.mod \
                error_handling.mod fitpack_utils.mod kinds.mod \
                nlps.mod qspl.mod system.mod zeroing_utils.mod

radin_utils.mod.f90:$(SRCDIR)/radin_utils.mod.F90
radin_utils.mod.o:radin_utils.mod.f90 error_handling.mod kinds.mod \
                parac.mod

ragg.mod.f90:   $(SRCDIR)/ragg.mod.F90
ragg.mod.o:     ragg.mod.f90 kinds.mod system.mod

raman_p_utils.mod.f90:$(SRCDIR)/raman_p_utils.mod.F90
raman_p_utils.mod.o:raman_p_utils.mod.f90 adat.mod cnst.mod \
                d_mat_p_utils.mod ddip.mod ddipo_utils.mod \
                dipomod.mod dotp_utils.mod elct.mod error_handling.mod \
                fft_maxfft.mod fileopen_utils.mod fileopenmod.mod \
                geq0mod.mod gvec.mod ions.mod kinds.mod kpnt.mod \
                mp_interface.mod nlps.mod opeigr_p_utils.mod \
                parac.mod response_pmod.mod rnlsm_utils.mod \
                rwfopt_p_utils.mod sfac.mod summat_utils.mod \
                system.mod timer.mod utils.mod zeroing_utils.mod

ranc_utils.mod.f90:$(SRCDIR)/ranc_utils.mod.F90
ranc_utils.mod.o:ranc_utils.mod.f90 ions.mod kinds.mod metr.mod \
                mp_interface.mod parac.mod prng_utils.mod setsc_utils.mod \
                system.mod

randtowf_utils.mod.f90:$(SRCDIR)/randtowf_utils.mod.F90
randtowf_utils.mod.o:randtowf_utils.mod.f90 cppt.mod geq0mod.mod \
                kinds.mod kpclean_utils.mod kpnt.mod kpts.mod \
                prng_utils.mod system.mod timer.mod

ranp_utils.mod.f90:$(SRCDIR)/ranp_utils.mod.F90
ranp_utils.mod.o:ranp_utils.mod.f90 cotr.mod ions.mod kinds.mod \
                parac.mod prng_utils.mod system.mod

ratom_utils.mod.f90:$(SRCDIR)/ratom_utils.mod.F90
ratom_utils.mod.o:ratom_utils.mod.f90 atom.mod atwf.mod clas.mod \
                cnst_dyn.mod coninp_utils.mod coor.mod cotr.mod \
                dpot.mod elct.mod error_handling.mod ghermit_utils.mod \
                inscan_utils.mod ions.mod kinds.mod linres.mod \
                meta_cell_utils.mod meta_colvar_inp_utils.mod \
                meta_dyn_def_utils.mod mm_dimmod.mod mm_input.mod \
                movi.mod mp_interface.mod mw.mod nlcc.mod nlps.mod \
                parac.mod pslo.mod ragg.mod readsr_utils.mod \
                recpnew_utils.mod recpupf_utils.mod rmas.mod \
                symm.mod system.mod tst2min_inp_utils.mod velocitinp_utils.mod \
                zeroing_utils.mod

rattle_utils.mod.f90:$(SRCDIR)/rattle_utils.mod.F90
rattle_utils.mod.o:rattle_utils.mod.f90 cnstfc_utils.mod cotr.mod \
                error_handling.mod kinds.mod mm_dim_utils.mod \
                mm_dimmod.mod parac.mod puttau_utils.mod system.mod \
                tpar.mod zeroing_utils.mod

rbfgs_utils.mod.f90:$(SRCDIR)/rbfgs_utils.mod.F90
rbfgs_utils.mod.o:rbfgs_utils.mod.f90 cotr.mod error_handling.mod \
                fixcom_utils.mod kinds.mod timer.mod

readff_utils.mod.f90:$(SRCDIR)/readff_utils.mod.F90
readff_utils.mod.o:readff_utils.mod.f90 clas.mod error_handling.mod \
                kinds.mod parac.mod readsr_utils.mod zeroing_utils.mod

readmod.mod.f90:$(SRCDIR)/readmod.mod.F90
readmod.mod.o:  readmod.mod.f90

read_prop_utils.mod.f90:$(SRCDIR)/read_prop_utils.mod.F90
read_prop_utils.mod.o:read_prop_utils.mod.f90 efld.mod error_handling.mod \
                inscan_utils.mod mp_interface.mod parac.mod \
                system.mod td_input.mod

readsr_utils.mod.f90:$(SRCDIR)/readsr_utils.mod.F90
readsr_utils.mod.o:readsr_utils.mod.f90 kinds.mod

readvan_utils.mod.f90:$(SRCDIR)/readvan_utils.mod.F90
readvan_utils.mod.o:readvan_utils.mod.f90 error_handling.mod \
                ions.mod kinds.mod parac.mod pslo.mod readsr_utils.mod \
                system.mod vdbp.mod vdbt.mod zeroing_utils.mod

recpnew_utils.mod.f90:$(SRCDIR)/recpnew_utils.mod.F90
recpnew_utils.mod.o:recpnew_utils.mod.f90 adat.mod atom.mod \
                cp_xc_utils.mod dpot.mod error_handling.mod \
                fitpack_utils.mod func.mod inscan_utils.mod \
                ions.mod kinds.mod machine.mod nlcc.mod nlps.mod \
                parac.mod pslo.mod ragg.mod readsr_utils.mod \
                readvan_utils.mod rmas.mod sgpp.mod special_functions.mod \
                system.mod vdbp.mod vdbt.mod

recpupf_utils.mod.f90:$(SRCDIR)/recpupf_utils.mod.F90
recpupf_utils.mod.o:recpupf_utils.mod.f90 adat.mod atom.mod \
                dpot.mod error_handling.mod inscan_utils.mod \
                ions.mod kinds.mod nlcc.mod nlps.mod pslo.mod \
                ragg.mod readsr_utils.mod recpnew_utils.mod \
                rmas.mod sgpp.mod system.mod vdbt.mod

reigs_utils.mod.f90:$(SRCDIR)/reigs_utils.mod.F90
reigs_utils.mod.o:reigs_utils.mod.f90 cnst.mod error_handling.mod \
                kinds.mod parac.mod ropt.mod spin.mod system.mod \
                timer.mod utils.mod

rekine_utils.mod.f90:$(SRCDIR)/rekine_utils.mod.F90
rekine_utils.mod.o:rekine_utils.mod.f90 dotp_utils.mod geq0mod.mod \
                harm.mod kinds.mod mp_interface.mod parac.mod \
                system.mod

repgen_utils.mod.f90:$(SRCDIR)/repgen_utils.mod.F90
repgen_utils.mod.o:repgen_utils.mod.f90 adat.mod cnst.mod error_handling.mod \
                fileopen_utils.mod fileopenmod.mod ions.mod \
                isos.mod kinds.mod movi.mod mp_interface.mod \
                parac.mod pimd.mod pinmtrans_utils.mod prng_utils.mod \
                store_types.mod system.mod

resetac_utils.mod.f90:$(SRCDIR)/resetac_utils.mod.F90
resetac_utils.mod.o:resetac_utils.mod.f90 kinds.mod system.mod \
                zeroing_utils.mod

reshaper.mod.f90:$(SRCDIR)/reshaper.mod.F90
reshaper.mod.o: reshaper.mod.f90 kinds.mod

respin_p_utils.mod.f90:$(SRCDIR)/respin_p_utils.mod.F90
respin_p_utils.mod.o:respin_p_utils.mod.f90 elct.mod error_handling.mod \
                inscan_utils.mod kinds.mod kpnt.mod kpts.mod \
                mp_interface.mod parac.mod readsr_utils.mod \
                response_pmod.mod store_types.mod system.mod

response_p.mod.f90:$(SRCDIR)/response_p.mod.F90
response_p.mod.o:response_p.mod.f90 kinds.mod

response_p_utils.mod.f90:$(SRCDIR)/response_p_utils.mod.F90
response_p_utils.mod.o:response_p_utils.mod.f90 coor.mod do_perturbation_p_utils.mod \
                elct.mod error_handling.mod filnmod.mod fnlalloc_utils.mod \
                kinds.mod kpts.mod mp_interface.mod parac.mod \
                response_pmod.mod ropt.mod rv30_utils.mod setirec_utils.mod \
                soft.mod store_types.mod system.mod testex_utils.mod \
                zeroing_utils.mod

restart_p_utils.mod.f90:$(SRCDIR)/restart_p_utils.mod.F90
restart_p_utils.mod.o:restart_p_utils.mod.f90 coor.mod error_handling.mod \
                fileopen_utils.mod fileopenmod.mod kinds.mod \
                mp_interface.mod parac.mod readsr_utils.mod \
                response_pmod.mod system.mod timer.mod zeroing_utils.mod

rgdiis_utils.mod.f90:$(SRCDIR)/rgdiis_utils.mod.F90
rgdiis_utils.mod.o:rgdiis_utils.mod.f90 cotr.mod error_handling.mod \
                fixcom_utils.mod kinds.mod odiis_utils.mod \
                system.mod timer.mod zeroing_utils.mod

rggen_utils.mod.f90:$(SRCDIR)/rggen_utils.mod.F90
rggen_utils.mod.o:rggen_utils.mod.f90 broy.mod cell.mod cppt.mod \
                error_handling.mod fileopen_utils.mod fileopenmod.mod \
                glopar_utils.mod gvec.mod kinds.mod latgen_utils.mod \
                loadpa_utils.mod metr.mod mp_interface.mod \
                parac.mod prmem_utils.mod setsc_utils.mod sphe.mod \
                system.mod timer.mod utils.mod

rgmopt_utils.mod.f90:$(SRCDIR)/rgmopt_utils.mod.F90
rgmopt_utils.mod.o:rgmopt_utils.mod.f90 zeroing_utils.mod kinds.mod \
                error_handling.mod timer.mod machine.mod mp_interface.mod \
                setsc_utils.mod dynit_utils.mod symtrz_utils.mod \
                chksym_utils.mod lsforce_utils.mod mm_qmmm_forcedr_utils.mod \
                wrccfl_utils.mod sdcell_utils.mod puttau_utils.mod \
                rrfo_utils.mod hessup_utils.mod hessout_utils.mod \
                hessin_utils.mod detdof_utils.mod sdion_utils.mod \
                adapttol_utils.mod cnstfc_utils.mod rgdiis_utils.mod \
                cnstpr_utils.mod rbfgs_utils.mod inr_dr_utils.mod \
                rlbfgs_utils.mod parac.mod atwf.mod spin.mod \
                ener.mod cnst.mod elct.mod tpar.mod pslo.mod \
                ions.mod soft.mod norm.mod ropt.mod coor.mod \
                sfac.mod cotr.mod nlcc.mod andr.mod andp.mod \
                nlps.mod fint.mod poin.mod kpts.mod kpnt.mod \
                store_types.mod metr.mod symm.mod xinr.mod \
                response_pmod.mod implhv.mod linres.mod bsym.mod \
                bsympnt.mod mm_dimmod.mod mm_input.mod lscal.mod \
                efld.mod cdftmod.mod system.mod vdwcmod.mod \
                hfxmod.mod wann.mod rinitwf_driver.mod setbsstate_utils.mod \
                mm_dim_utils.mod bs_forces_diag_utils.mod newcell_utils.mod \
                totstr_utils.mod dum2_utils.mod moverho_utils.mod \
                copot_utils.mod updrho_utils.mod forces_diag_utils.mod \
                calc_alm_utils.mod extrap_utils.mod localize_utils.mod \
                lr_tddft_utils.mod updwf_utils.mod cdft_utils.mod \
                k_updwf_utils.mod ortho_utils.mod rhopri_utils.mod \
                finalp_utils.mod phfac_utils.mod wrener_utils.mod \
                forcep_utils.mod wrener_utils.mod wrgeo_utils.mod \
                gsize_utils.mod hesele_utils.mod geofile_utils.mod \
                setirec_utils.mod wv30_utils.mod initrun_driver.mod \
                initrun_utils.mod testex_utils.mod rnlsm_utils.mod \
                forcedr_driver.mod zeroing_utils.mod kinds.mod \
                error_handling.mod timer.mod newcell_utils.mod \
                parac.mod elct.mod nlcc.mod store_types.mod \
                pslo.mod fint.mod atwf.mod system.mod vdwcmod.mod \
                rgdiis_utils.mod moverho_utils.mod rbfgs_utils.mod \
                rlbfgs_utils.mod rrfo_utils.mod sdion_utils.mod \
                copot_utils.mod ddipo_utils.mod forces_diag_utils.mod \
                calc_alm_utils.mod inr_dr_utils.mod lr_tddft_utils.mod \
                updwf_utils.mod ortho_utils.mod rhopri_utils.mod \
                rhoofr_utils.mod initrun_utils.mod rnlsm_utils.mod \
                forcedr_utils.mod

rgs_utils.mod.f90:$(SRCDIR)/rgs_utils.mod.F90
rgs_utils.mod.o:rgs_utils.mod.f90 cp_grp_utils.mod error_handling.mod \
                geq0mod.mod kinds.mod mp_interface.mod parac.mod \
                state_utils.mod system.mod timer.mod utils.mod \
                zeroing_utils.mod

rgsvan_utils.mod.f90:$(SRCDIR)/rgsvan_utils.mod.F90
rgsvan_utils.mod.o:rgsvan_utils.mod.f90 csmat_utils.mod kinds.mod \
                rgs_utils.mod sfac.mod system.mod timer.mod

rho1ofr_utils.mod.f90:$(SRCDIR)/rho1ofr_utils.mod.F90
rho1ofr_utils.mod.o:rho1ofr_utils.mod.f90 cnst.mod cppt.mod \
                dotp_utils.mod ener.mod error_handling.mod \
                fft.mod fftmain_utils.mod geq0mod.mod kinds.mod \
                mp_interface.mod parac.mod pslo.mod reshaper.mod \
                spin.mod system.mod timer.mod zeroing_utils.mod

rho1pri_utils.mod.f90:$(SRCDIR)/rho1pri_utils.mod.F90
rho1pri_utils.mod.o:rho1pri_utils.mod.f90 cppt.mod elct.mod \
                error_handling.mod fftmain_utils.mod forcep_utils.mod \
                kinds.mod parac.mod rho1ofr_utils.mod system.mod

rhodiis_utils.mod.f90:$(SRCDIR)/rhodiis_utils.mod.F90
rhodiis_utils.mod.o:rhodiis_utils.mod.f90 andr.mod kinds.mod \
                mp_interface.mod odiis_utils.mod parac.mod \
                zeroing_utils.mod

rhoofr_c_utils.mod.f90:$(SRCDIR)/rhoofr_c_utils.mod.F90
rhoofr_c_utils.mod.o:rhoofr_c_utils.mod.f90 augchg_utils.mod \
                cp_grp_utils.mod density_utils.mod elct.mod \
                ener.mod error_handling.mod fftmain_utils.mod \
                fftnew_utils.mod ions.mod kinds.mod kpnt.mod \
                kpts.mod moverho_utils.mod mp_interface.mod \
                parac.mod part_1d.mod prcp.mod pslo.mod reshaper.mod \
                rhov_utils.mod ropt.mod rswfmod.mod sfac.mod \
                special_functions.mod spin.mod state_utils.mod \
                symm.mod symtrz_utils.mod system.mod timer.mod \
                zeroing_utils.mod

rhoofr_kdp_utils.mod.f90:$(SRCDIR)/rhoofr_kdp_utils.mod.F90
rhoofr_kdp_utils.mod.o:rhoofr_kdp_utils.mod.f90 augchg_utils.mod \
                cnst.mod cppt.mod dotp_utils.mod elct.mod ener.mod \
                error_handling.mod fft.mod fft_maxfft.mod fftmain_utils.mod \
                geq0mod.mod ions.mod kinds.mod moverho_utils.mod \
                mp_interface.mod parac.mod pslo.mod reshaper.mod \
                rhov_utils.mod ropt.mod sfac.mod spin.mod symm.mod \
                symtrz_utils.mod system.mod timer.mod zeroing_utils.mod

rhoofr_p_utils.mod.f90:$(SRCDIR)/rhoofr_p_utils.mod.F90
rhoofr_p_utils.mod.o:rhoofr_p_utils.mod.f90 cnst.mod cppt.mod \
                dotp_utils.mod elct.mod ener.mod error_handling.mod \
                fft.mod fft_maxfft.mod fftmain_utils.mod geq0mod.mod \
                ions.mod kinds.mod moverho_utils.mod mp_interface.mod \
                parac.mod perturbation_p_utils.mod pslo.mod \
                reshaper.mod response_pmod.mod rhov_utils.mod \
                ropt.mod spin.mod symtrz_utils.mod system.mod \
                timer.mod zeroing_utils.mod

rhoofr_utils.mod.f90:$(SRCDIR)/rhoofr_utils.mod.F90
rhoofr_utils.mod.o:rhoofr_utils.mod.f90 augchg_utils.mod cnst.mod \
                cp_cuda_types.mod cp_cudensity_utils.mod cp_cufft_types.mod \
                cp_curho_types.mod cp_curho_utils.mod cp_cuwfn_types.mod \
                cp_grp_utils.mod cppt.mod cuda_types.mod cuda_utils.mod \
                cuuser_utils.mod density_utils.mod dg.mod elct.mod \
                ener.mod error_handling.mod fft.mod fftmain_utils.mod \
                fftnew_utils.mod geq0mod.mod ions.mod kin_energy_utils.mod \
                kinds.mod moverho_utils.mod mp_interface.mod \
                parac.mod part_1d.mod pslo.mod rho1ofr_utils.mod \
                rhov_utils.mod ropt.mod rswfmod.mod sfac.mod \
                spin.mod state_utils.mod symm.mod symtrz_utils.mod \
                system.mod tauofr_utils.mod thread_view_types.mod \
                thread_view_utils.mod timer.mod zeroing_utils.mod \
                nvtx_utils.mod

rhopri_utils.mod.f90:$(SRCDIR)/rhopri_utils.mod.F90
rhopri_utils.mod.o:rhopri_utils.mod.f90 bsym.mod cdftmod.mod \
                cppt.mod eicalc_utils.mod elf_utils.mod elstpo_utils.mod \
                ener.mod error_handling.mod fftmain_utils.mod \
                hip_utils.mod kinds.mod kpts.mod lsd_elf_utils.mod \
                meta_multiple_walkers_utils.mod mw.mod noforce_utils.mod \
                norhoe_utils.mod parac.mod phfac_utils.mod \
                pimd.mod prden.mod readsr_utils.mod response_pmod.mod \
                rhoofr_c_utils.mod rhoofr_utils.mod rnlsm_utils.mod \
                ropt.mod spin.mod store_types.mod system.mod \
                timer.mod utils.mod zeroing_utils.mod

rhov1_utils.mod.f90:$(SRCDIR)/rhov1_utils.mod.F90
rhov1_utils.mod.o:rhov1_utils.mod.f90 cppt.mod error_handling.mod \
                ions.mod kinds.mod nlps.mod pslo.mod sfac.mod \
                system.mod timer.mod zeroing_utils.mod

rhov_utils.mod.f90:$(SRCDIR)/rhov_utils.mod.F90
rhov_utils.mod.o:rhov_utils.mod.f90 cppt.mod elct.mod error_handling.mod \
                fftmain_utils.mod fftnew_utils.mod geq0mod.mod \
                ions.mod kinds.mod nlps.mod parac.mod pslo.mod \
                qvan2_utils.mod sfac.mod system.mod timer.mod \
                zeroing_utils.mod

rinforce_nuc_utils.mod.f90:$(SRCDIR)/rinforce_nuc_utils.mod.F90
rinforce_nuc_utils.mod.o:rinforce_nuc_utils.mod.f90 aainit_utils.mod \
                aavan.mod atom.mod cnst.mod cppt.mod cvan.mod \
                dpot.mod eam.mod eam_pot_utils.mod ehpsi_utils.mod \
                ener.mod error_handling.mod fint.mod formf_utils.mod \
                ions.mod kinds.mod kpts.mod mm_input.mod mm_ion_dens.mod \
                nlcc.mod nlccset_utils.mod nlps.mod parac.mod \
                prmem_utils.mod pslo.mod qrada_s_utils.mod \
                qspl.mod ragg.mod rinforce_utils.mod rnlin_utils.mod \
                rnlset_utils.mod sgpp.mod str2.mod system.mod \
                timer.mod vdbinit_utils.mod vdbp.mod zeroing_utils.mod

rinforce_utils.mod.f90:$(SRCDIR)/rinforce_utils.mod.F90
rinforce_utils.mod.o:rinforce_utils.mod.f90 aainit_utils.mod \
                aavan.mod atom.mod cnst.mod cppt.mod cvan.mod \
                dpot.mod eam.mod eam_pot_utils.mod ehpsi_utils.mod \
                elct.mod ener.mod error_handling.mod fint.mod \
                fitpack_utils.mod formf_utils.mod geq0mod.mod \
                gvec.mod ions.mod isos.mod kinds.mod kpclean_utils.mod \
                kpnt.mod kpts.mod mm_input.mod mm_ion_dens.mod \
                mp_interface.mod nlcc.mod nlccset_utils.mod \
                nlps.mod parac.mod prmem_utils.mod pslo.mod \
                qrada_s_utils.mod qspl.mod ragg.mod rnlin_utils.mod \
                rnlset_utils.mod sgpp.mod sphe.mod str2.mod \
                system.mod timer.mod utils.mod vdbinit_utils.mod \
                vdbp.mod ylmr2_utils.mod zeroing_utils.mod

rinit_utils.mod.f90:$(SRCDIR)/rinit_utils.mod.F90
rinit_utils.mod.o:rinit_utils.mod.f90 broy.mod cell.mod clas.mod \
                cnst.mod elct.mod error_handling.mod fint.mod \
                gvec.mod isos.mod kdpc.mod kdpoints_utils.mod \
                kinds.mod kpnt.mod kpts.mod metr.mod nlps.mod \
                parac.mod prcp.mod readsr_utils.mod response_pmod.mod \
                rggen_utils.mod rkpnt_utils.mod sphe.mod store_types.mod \
                symm.mod system.mod timer.mod zeroing_utils.mod

rinitwf_driver.mod.f90:$(SRCDIR)/rinitwf_driver.mod.F90
rinitwf_driver.mod.o:rinitwf_driver.mod.f90 atomwf_utils.mod \
                atwf.mod copot_utils.mod dotp_utils.mod error_handling.mod \
                forcedr_driver.mod ions.mod kinds.mod kpts.mod \
                mm_input.mod mp_interface.mod newcell_utils.mod \
                nlcc.mod ortho_utils.mod parac.mod phfac_utils.mod \
                pslo.mod randtowf_utils.mod reshaper.mod rnlsm_utils.mod \
                ropt.mod setbasis_utils.mod sphe.mod system.mod \
                timer.mod tpar.mod utils.mod

rinitwf_utils.mod.f90:$(SRCDIR)/rinitwf_utils.mod.F90
rinitwf_utils.mod.o:rinitwf_utils.mod.f90 atomwf_utils.mod \
                copot_utils.mod forcedr_utils.mod newcell_utils.mod \
                ortho_utils.mod system.mod timer.mod

rinvel_utils.mod.f90:$(SRCDIR)/rinvel_utils.mod.F90
rinvel_utils.mod.o:rinvel_utils.mod.f90 cnst.mod cnst_dyn.mod \
                coor.mod cotr.mod ekinpp_utils.mod ions.mod \
                kinds.mod metr.mod mm_input.mod mp_interface.mod \
                nose.mod parac.mod pimd.mod prng_utils.mod \
                puttau_utils.mod rmas.mod system.mod zeroing_utils.mod

rk4ov_utils.mod.f90:$(SRCDIR)/rk4ov_utils.mod.F90
rk4ov_utils.mod.o:rk4ov_utils.mod.f90 kinds.mod shop.mod shop_rest.mod

rkpnt_utils.mod.f90:$(SRCDIR)/rkpnt_utils.mod.F90
rkpnt_utils.mod.o:rkpnt_utils.mod.f90 coor.mod cppt.mod envj.mod \
                error_handling.mod fileopen_utils.mod fileopenmod.mod \
                gvec.mod ions.mod k290_utils.mod kinds.mod \
                kpnt.mod kpts.mod machine.mod mp_interface.mod \
                parac.mod prmem_utils.mod prng_utils.mod readsr_utils.mod \
                sphe.mod symm.mod system.mod timer.mod zeroing_utils.mod

rlbfgs_io.mod.f90:$(SRCDIR)/rlbfgs_io.mod.F90
rlbfgs_io.mod.o:rlbfgs_io.mod.f90 cotr.mod error_handling.mod \
                lscal.mod parac.mod utils.mod

rlbfgs_utils.mod.f90:$(SRCDIR)/rlbfgs_utils.mod.F90
rlbfgs_utils.mod.o:rlbfgs_utils.mod.f90 adapttol_utils.mod \
                coor.mod cotr.mod error_handling.mod hessin_utils.mod \
                ions.mod kinds.mod lscal.mod parac.mod sdion_utils.mod \
                secder_utils.mod store_types.mod system.mod \
                timer.mod vibana_utils.mod zeroing_utils.mod \
                kinds.mod error_handling.mod timer.mod mp_interface.mod \
                parac.mod lscal.mod rlbfgs_utils.mod

rmas.mod.f90:   $(SRCDIR)/rmas.mod.F90
rmas.mod.o:     rmas.mod.f90 kinds.mod system.mod

rnl_dk_p_utils.mod.f90:$(SRCDIR)/rnl_dk_p_utils.mod.F90
rnl_dk_p_utils.mod.o:rnl_dk_p_utils.mod.f90 error_handling.mod \
                geq0mod.mod ions.mod kinds.mod kpnt.mod nlps.mod \
                response_pmod.mod sumfnl_utils.mod system.mod \
                timer.mod zeroing_utils.mod

rnlfl_utils.mod.f90:$(SRCDIR)/rnlfl_utils.mod.F90
rnlfl_utils.mod.o:rnlfl_utils.mod.f90 cvan.mod error_handling.mod \
                ions.mod kinds.mod nlps.mod parac.mod pslo.mod \
                sfac.mod system.mod timer.mod

rnlfor_utils.mod.f90:$(SRCDIR)/rnlfor_utils.mod.F90
rnlfor_utils.mod.o:rnlfor_utils.mod.f90 cvan.mod error_handling.mod \
                ions.mod kinds.mod mp_interface.mod nlps.mod \
                parac.mod pslo.mod sfac.mod sgpp.mod spin.mod \
                system.mod timer.mod zeroing_utils.mod

rnlin_utils.mod.f90:$(SRCDIR)/rnlin_utils.mod.F90
rnlin_utils.mod.o:rnlin_utils.mod.f90 atom.mod bessm_utils.mod \
                cnst.mod dcacp_utils.mod dpot.mod error_handling.mod \
                fitpack_utils.mod kinds.mod mp_interface.mod \
                nlps.mod parac.mod pslo.mod qspl.mod sgpp.mod \
                system.mod timer.mod utils.mod ylmr_utils.mod \
                zeroing_utils.mod

rnlrh_utils.mod.f90:$(SRCDIR)/rnlrh_utils.mod.F90
rnlrh_utils.mod.o:rnlrh_utils.mod.f90 cvan.mod elct.mod ener.mod \
                error_handling.mod ions.mod kinds.mod kpnt.mod \
                nlps.mod parac.mod pslo.mod sfac.mod sgpp.mod \
                spin.mod system.mod timer.mod

rnlset_utils.mod.f90:$(SRCDIR)/rnlset_utils.mod.F90
rnlset_utils.mod.o:rnlset_utils.mod.f90 atom.mod dcacp_utils.mod \
                dpot.mod error_handling.mod fitpack_utils.mod \
                kinds.mod nlps.mod sgpp.mod system.mod timer.mod \
                utils.mod

rnlsm1_utils.mod.f90:$(SRCDIR)/rnlsm1_utils.mod.F90
rnlsm1_utils.mod.o:rnlsm1_utils.mod.f90 cppt.mod error_handling.mod \
                geq0mod.mod ions.mod kinds.mod kpnt.mod kpts.mod \
                mm_dimmod.mod mp_interface.mod nlps.mod nvtx_utils.mod \
                parac.mod part_1d.mod sfac.mod sumfnl_utils.mod \
                system.mod timer.mod zeroing_utils.mod

rnlsm_2d_utils.mod.f90:$(SRCDIR)/rnlsm_2d_utils.mod.F90
rnlsm_2d_utils.mod.o:rnlsm_2d_utils.mod.f90 cppt.mod error_handling.mod \
                geq0mod.mod ions.mod kinds.mod mp_interface.mod \
                nlps.mod parac.mod sfac.mod system.mod timer.mod

rnlsm2_utils.mod.f90:$(SRCDIR)/rnlsm2_utils.mod.F90
rnlsm2_utils.mod.o:rnlsm2_utils.mod.f90 cppt.mod error_handling.mod \
                geq0mod.mod ions.mod kinds.mod kpnt.mod kpts.mod \
                mm_dimmod.mod mp_interface.mod nlps.mod parac.mod \
                sfac.mod system.mod timer.mod zeroing_utils.mod

rnlsmd_utils.mod.f90:$(SRCDIR)/rnlsmd_utils.mod.F90
rnlsmd_utils.mod.o:rnlsmd_utils.mod.f90 cppt.mod error_handling.mod \
                geq0mod.mod ions.mod kinds.mod kpts.mod mm_dimmod.mod \
                nlps.mod parac.mod sfac.mod sumfnl_utils.mod \
                system.mod timer.mod zeroing_utils.mod

rnlsm_p_utils.mod.f90:$(SRCDIR)/rnlsm_p_utils.mod.F90
rnlsm_p_utils.mod.o:rnlsm_p_utils.mod.f90 cppt.mod error_handling.mod \
                geq0mod.mod ions.mod kinds.mod kpts.mod mp_interface.mod \
                nlps.mod parac.mod sfac.mod system.mod timer.mod \
                zeroing_utils.mod

rnlsm_utils.mod.f90:$(SRCDIR)/rnlsm_utils.mod.F90
rnlsm_utils.mod.o:rnlsm_utils.mod.f90 kinds.mod nlps.mod rnlsm1_utils.mod \
                rnlsm2_utils.mod rnlsmd_utils.mod system.mod \
                timer.mod

ropt.mod.f90:   $(SRCDIR)/ropt.mod.F90
ropt.mod.o:     ropt.mod.f90

rortog_utils.mod.f90:$(SRCDIR)/rortog_utils.mod.F90
rortog_utils.mod.o:rortog_utils.mod.f90 elct.mod error_handling.mod \
                jrotation_utils.mod kinds.mod mp_interface.mod \
                ovlap_utils.mod parac.mod reigs_utils.mod spin.mod \
                system.mod timer.mod tpar.mod utils.mod zeroing_utils.mod \
                kinds.mod error_handling.mod timer.mod mp_interface.mod \
                system.mod parac.mod spin.mod ovlap_utils.mod \
                summat_utils.mod utils.mod rotate_utils.mod

rortv_utils.mod.f90:$(SRCDIR)/rortv_utils.mod.F90
rortv_utils.mod.o:rortv_utils.mod.f90 crotwf_utils.mod dotp_utils.mod \
                error_handling.mod geq0mod.mod harm.mod jrotation_utils.mod \
                kinds.mod linalg_utils.mod mp_interface.mod \
                nort.mod ovlap_utils.mod parac.mod rotate_utils.mod \
                spin.mod system.mod tpar.mod zeroing_utils.mod

rotate_my_wannier_manno_p_utils.mod.f90:$(SRCDIR)/rotate_my_wannier_manno_p_utils.mod.F90
rotate_my_wannier_manno_p_utils.mod.o:rotate_my_wannier_manno_p_utils.mod.f90 \
                coor.mod cppt.mod ener.mod error_handling.mod \
                fileopen_utils.mod fileopenmod.mod fnlalloc_utils.mod \
                forcep_utils.mod forces_driver.mod forces_utils.mod \
                gndstate_p_utils.mod kinds.mod legendre_p_utils.mod \
                metr.mod ovlap_utils.mod parac.mod phfac_utils.mod \
                readsr_utils.mod response_pmod.mod symm.mod \
                system.mod timer.mod utils.mod zeroing_utils.mod

rotate_my_wannier_para_p_utils.mod.f90:$(SRCDIR)/rotate_my_wannier_para_p_utils.mod.F90
rotate_my_wannier_para_p_utils.mod.o:rotate_my_wannier_para_p_utils.mod.f90 \
                coor.mod cppt.mod ener.mod error_handling.mod \
                fileopen_utils.mod fileopenmod.mod fnlalloc_utils.mod \
                forcep_utils.mod forces_driver.mod forces_utils.mod \
                geq0mod.mod gndstate_p_utils.mod kinds.mod \
                lowdin_utils.mod metr.mod mp_interface.mod \
                ovlap_utils.mod parac.mod phfac_utils.mod readsr_utils.mod \
                response_pmod.mod symm.mod system.mod timer.mod \
                utils.mod zeroing_utils.mod

rotate_utils.mod.f90:$(SRCDIR)/rotate_utils.mod.F90
rotate_utils.mod.o:rotate_utils.mod.f90 error_handling.mod \
                kinds.mod nvtx_utils.mod timer.mod kinds.mod \
                error_handling.mod timer.mod mp_interface.mod \
                cuda_types.mod cuda_utils.mod cublas_types.mod \
                cublas_utils.mod sizeof_kinds.mod cp_cuwfn_types.mod \
                nvtx_utils.mod kinds.mod error_handling.mod \
                timer.mod system.mod parac.mod spin.mod

rotvel_utils.mod.f90:$(SRCDIR)/rotvel_utils.mod.F90
rotvel_utils.mod.o:rotvel_utils.mod.f90 cnst.mod ekinpp_utils.mod \
                ions.mod kinds.mod nose.mod rmas.mod system.mod \
                zeroing_utils.mod

rpiiint_utils.mod.f90:$(SRCDIR)/rpiiint_utils.mod.F90
rpiiint_utils.mod.o:rpiiint_utils.mod.f90 cnst.mod eam.mod \
                eam_pot_utils.mod ions.mod kinds.mod metr.mod \
                mp_interface.mod parac.mod pbc_utils.mod ragg.mod \
                special_functions.mod system.mod timer.mod

rrandd_utils.mod.f90:$(SRCDIR)/rrandd_utils.mod.F90
rrandd_utils.mod.o:rrandd_utils.mod.f90 kinds.mod parac.mod \
                prng_utils.mod spin.mod system.mod

rrane_utils.mod.f90:$(SRCDIR)/rrane_utils.mod.F90
rrane_utils.mod.o:rrane_utils.mod.f90 cppt.mod error_handling.mod \
                geq0mod.mod gvec.mod kinds.mod kpclean_utils.mod \
                kpts.mod prng_utils.mod pslo.mod system.mod \
                utils.mod

rreadf_utils.mod.f90:$(SRCDIR)/rreadf_utils.mod.F90
rreadf_utils.mod.o:rreadf_utils.mod.f90 cnst.mod error_handling.mod \
                fileopen_utils.mod fileopenmod.mod ions.mod \
                kinds.mod mp_interface.mod parac.mod pimd.mod \
                readsr_utils.mod system.mod velupi_utils.mod

rrfo_utils.mod.f90:$(SRCDIR)/rrfo_utils.mod.F90
rrfo_utils.mod.o:rrfo_utils.mod.f90 cotr.mod error_handling.mod \
                fixcom_utils.mod kinds.mod parac.mod system.mod \
                timer.mod utils.mod zeroing_utils.mod

rscpot_utils.mod.f90:$(SRCDIR)/rscpot_utils.mod.F90
rscpot_utils.mod.o:rscpot_utils.mod.f90 cdft_utils.mod cnst_dyn.mod \
                elct.mod ener.mod epr_efg_utils.mod hubbardu.mod \
                hubbardu_utils.mod kinds.mod kpnt.mod kpts.mod \
                mm_input.mod mp_interface.mod parac.mod prop.mod \
                rhoofr_c_utils.mod rhoofr_utils.mod rnlfor_utils.mod \
                rnlrh_utils.mod ropt.mod rpiiint_utils.mod \
                spin.mod stress_utils.mod system.mod timer.mod \
                vdw_utils.mod vdwcmod.mod vofrho_utils.mod \
                vofrhoc_utils.mod zeroing_utils.mod

rscve_utils.mod.f90:$(SRCDIR)/rscve_utils.mod.F90
rscve_utils.mod.o:rscve_utils.mod.f90 geq0mod.mod kinds.mod \
                parac.mod

rscvp_utils.mod.f90:$(SRCDIR)/rscvp_utils.mod.F90
rscvp_utils.mod.o:rscvp_utils.mod.f90 cnst_dyn.mod ions.mod \
                kinds.mod parac.mod system.mod

rswf.mod.f90:   $(SRCDIR)/rswf.mod.F90
rswf.mod.o:     rswf.mod.f90 kinds.mod

rv30_utils.mod.f90:$(SRCDIR)/rv30_utils.mod.F90
rv30_utils.mod.o:rv30_utils.mod.f90 cdftmod.mod cell.mod clas.mod \
                cotr.mod elct.mod ener.mod error_handling.mod \
                fileopenmod.mod filnmod.mod geq0mod.mod glemod.mod \
                io_utils.mod ions.mod kinds.mod kpnt.mod kpts.mod \
                lscal.mod machine.mod meta_multiple_walkers_utils.mod \
                metr.mod mm_dimmod.mod mm_extrap.mod mp_interface.mod \
                mw.mod nose.mod parac.mod pimd.mod poin.mod \
                prng.mod prng_utils.mod readsr_utils.mod rlbfgs_io.mod \
                rw_linres_utils.mod shop_rest.mod shop_rest_2.mod \
                spin.mod store_types.mod string_utils.mod symm.mod \
                system.mod timer.mod utils.mod wfnio_utils.mod \
                bicanonicalCpmd.mod zeroing_utils.mod kinds.mod \
                error_handling.mod timer.mod mp_interface.mod \
                system.mod parac.mod zeroing_utils.mod

rwfopt_nuc_utils.mod.f90:$(SRCDIR)/rwfopt_nuc_utils.mod.F90
rwfopt_nuc_utils.mod.o:rwfopt_nuc_utils.mod.f90 andp.mod coor.mod \
                dynit_utils.mod efld.mod elct.mod enbandpri_utils.mod \
                ener.mod error_handling.mod finalp_utils.mod \
                geofile_utils.mod gsize_utils.mod initrun_utils.mod \
                kinds.mod kpts.mod machine.mod mm_dim_utils.mod \
                mm_dimmod.mod mm_input.mod norm.mod parac.mod \
                poin.mod prmem_utils.mod rhopri_utils.mod ropt.mod \
                rscpot_utils.mod spin.mod store_types.mod system.mod \
                timer.mod totstr_utils.mod updrho_utils.mod \
                updwf_utils.mod wrener_utils.mod wv30_utils.mod \
                zeroing_utils.mod kinds.mod error_handling.mod \
                timer.mod fft_maxfft.mod system.mod parac.mod \
                spin.mod ener.mod cnst.mod rscpot_utils.mod

rwfopt_p_utils.mod.f90:$(SRCDIR)/rwfopt_p_utils.mod.F90
rwfopt_p_utils.mod.o:rwfopt_p_utils.mod.f90 coor.mod dotp_utils.mod \
                dynit_utils.mod elct.mod ener.mod error_handling.mod \
                fft_maxfft.mod fileopen_utils.mod fileopenmod.mod \
                hesele_p_utils.mod hpsi_utils.mod kinds.mod \
                machine.mod mp_interface.mod norm.mod ovlap_utils.mod \
                parac.mod perturbation_p_utils.mod response_pmod.mod \
                rho1ofr_utils.mod ropt.mod rpiiint_utils.mod \
                simple_model_p_utils.mod soft.mod spin.mod \
                system.mod testex_utils.mod tpar.mod updwf_p_utils.mod \
                utils.mod vofrho_utils.mod vpsi_p_utils.mod \
                wrener_utils.mod zeroing_utils.mod

rwfopt_utils.mod.f90:$(SRCDIR)/rwfopt_utils.mod.F90
rwfopt_utils.mod.o:rwfopt_utils.mod.f90 andp.mod cdft_utils.mod \
                cdftmod.mod cnst.mod coor.mod cplngs_utils.mod \
                cplngsmod.mod ddipo_utils.mod detdof_utils.mod \
                dynit_utils.mod efld.mod ehrenfest_utils.mod \
                elct.mod enbandpri_utils.mod ener.mod epr_efg_utils.mod \
                error_handling.mod fft_maxfft.mod fileopen_utils.mod \
                finalp_utils.mod forcep_utils.mod forces_diag_utils.mod \
                geofile_utils.mod gsize_utils.mod hfxmod.mod \
                hubbardu.mod hubbardu_utils.mod initrun_driver.mod \
                initrun_utils.mod ions.mod isos.mod k_updwf_utils.mod \
                kinds.mod kpts.mod linres.mod localize_utils.mod \
                locpot.mod lscal.mod machine.mod mm_dim_utils.mod \
                mm_dimmod.mod mm_input.mod mm_qmmm_forcedr_utils.mod \
                mp_interface.mod norm.mod parac.mod poin.mod \
                prmem_utils.mod rhoofr_c_utils.mod rhoofr_utils.mod \
                rhopri_utils.mod rinitwf_driver.mod ropt.mod \
                rscpot_utils.mod rswfmod.mod rv30_utils.mod \
                setirec_utils.mod spin.mod store_types.mod \
                syscomb_utils.mod system.mod td_input.mod td_utils.mod \
                timer.mod totstr_utils.mod transme_utils.mod \
                updrho_utils.mod updwf_utils.mod utils.mod \
                vdwcmod.mod wann.mod wrener_utils.mod wv30_utils.mod \
                zeroing_utils.mod

rw_linres_utils.mod.f90:$(SRCDIR)/rw_linres_utils.mod.F90
rw_linres_utils.mod.o:rw_linres_utils.mod.f90 error_handling.mod \
                io_utils.mod kinds.mod linres.mod mp_interface.mod \
                parac.mod system.mod utils.mod zeroing_utils.mod

rwswap_utils.mod.f90:$(SRCDIR)/rwswap_utils.mod.F90
rwswap_utils.mod.o:rwswap_utils.mod.f90 utils.mod zeroing_utils.mod \
                kinds.mod error_handling.mod timer.mod system.mod \
                parac.mod elct.mod cppt.mod nlps.mod fint.mod \
                swap.mod kpts.mod kpnt.mod sphe.mod phfac_utils.mod \
                kinds.mod error_handling.mod timer.mod rinforce_utils.mod \
                rkpnt_utils.mod ehpsi_utils.mod calc_alm_utils.mod \
                system.mod parac.mod cppt.mod nlps.mod fint.mod \
                kpts.mod kpnt.mod sphe.mod swap.mod phfac_utils.mod \
                kinds.mod error_handling.mod timer.mod system.mod \
                parac.mod fint.mod kpts.mod sphe.mod phfac_utils.mod \
                zeroing_utils.mod kinds.mod error_handling.mod \
                timer.mod readsr_utils.mod system.mod parac.mod \
                envj.mod kpts.mod swap.mod zeroing_utils.mod \
                kinds.mod error_handling.mod timer.mod swap.mod \
                parac.mod kinds.mod error_handling.mod timer.mod \
                swap.mod parac.mod utils.mod kinds.mod error_handling.mod \
                timer.mod swap.mod parac.mod kinds.mod error_handling.mod \
                timer.mod machine.mod swap.mod parac.mod kinds.mod \
                error_handling.mod timer.mod swap.mod parac.mod \
                zeroing_utils.mod kinds.mod error_handling.mod \
                timer.mod readsr_utils.mod swap.mod parac.mod \
                kinds.mod error_handling.mod timer.mod swap.mod \
                parac.mod kinds.mod error_handling.mod timer.mod \
                system.mod parac.mod kpts.mod swap.mod kinds.mod \
                error_handling.mod timer.mod swap.mod kpts.mod \
                kinds.mod error_handling.mod timer.mod readsr_utils.mod \
                swap.mod parac.mod kinds.mod error_handling.mod \
                timer.mod swap.mod parac.mod readsr_utils.mod

sample_utils.mod.f90:$(SRCDIR)/sample_utils.mod.F90
sample_utils.mod.o:sample_utils.mod.f90 error_handling.mod \
                fileopenmod.mod filnmod.mod machine.mod mp_interface.mod \
                parac.mod readsr_utils.mod store_types.mod \
                system.mod timer.mod

saop_utils.mod.f90:$(SRCDIR)/saop_utils.mod.F90
saop_utils.mod.o:saop_utils.mod.f90 cnst.mod error_handling.mod \
                functionals_utils.mod kinds.mod lsd_func_utils.mod \
                rho1ofr_utils.mod spin.mod system.mod zeroing_utils.mod

scrp.mod.f90:   $(SRCDIR)/scrp.mod.F90
scrp.mod.o:     scrp.mod.f90 kinds.mod

sdcell_utils.mod.f90:$(SRCDIR)/sdcell_utils.mod.F90
sdcell_utils.mod.o:sdcell_utils.mod.f90 ions.mod kinds.mod \
                metr.mod parac.mod prcp.mod setsc_utils.mod \
                system.mod tpar.mod

sd_ii_utils.mod.f90:$(SRCDIR)/sd_ii_utils.mod.F90
sd_ii_utils.mod.o:sd_ii_utils.mod.f90 cnst.mod ions.mod kinds.mod \
                metr.mod parac.mod pbc_utils.mod ragg.mod special_functions.mod \
                system.mod timer.mod

sdion_utils.mod.f90:$(SRCDIR)/sdion_utils.mod.F90
sdion_utils.mod.o:sdion_utils.mod.f90 cotr.mod ener.mod error_handling.mod \
                fixcom_utils.mod ions.mod kinds.mod parac.mod \
                sort_utils.mod tpar.mod zeroing_utils.mod

sdlinres_utils.mod.f90:$(SRCDIR)/sdlinres_utils.mod.F90
sdlinres_utils.mod.o:sdlinres_utils.mod.f90 adat.mod canon_utils.mod \
                coor.mod copot_utils.mod cotr.mod detdof_utils.mod \
                dynit_utils.mod eicalc_utils.mod eind_ii_utils.mod \
                eind_loc_utils.mod eind_nl_utils.mod elct.mod \
                ener.mod error_handling.mod fileopen_utils.mod \
                fileopenmod.mod fnlalloc_utils.mod forcedr_driver.mod \
                forcedr_utils.mod forcep_utils.mod hessin_utils.mod \
                hessout_utils.mod initrun_driver.mod initrun_utils.mod \
                ions.mod isos.mod kinds.mod kpts.mod ksdiag_utils.mod \
                linres.mod lr_in_utils.mod lr_xcpot_utils.mod \
                machine.mod mp_interface.mod nl_res_utils.mod \
                nlcc.mod nlps.mod norm.mod opt_lr_utils.mod \
                ortho_utils.mod parac.mod poin.mod rho1ofr_utils.mod \
                rhoofr_utils.mod rmas.mod rnlsm_2d_utils.mod \
                rnlsm_utils.mod ropt.mod sd_ii_utils.mod sd_loc2_utils.mod \
                sd_loc_utils.mod sd_nl2_utils.mod sd_nl_utils.mod \
                secder_utils.mod setirec_utils.mod sfac.mod \
                soft.mod spin.mod store_types.mod symm.mod \
                symtrz_utils.mod system.mod testex_utils.mod \
                updwf_utils.mod utils.mod wrener_utils.mod \
                wrgeo_utils.mod wv30_utils.mod zeroing_utils.mod

sd_loc2_utils.mod.f90:$(SRCDIR)/sd_loc2_utils.mod.F90
sd_loc2_utils.mod.o:sd_loc2_utils.mod.f90 cppt.mod fftmain_utils.mod \
                geq0mod.mod ions.mod kinds.mod parac.mod sfac.mod \
                system.mod

sd_loc_utils.mod.f90:$(SRCDIR)/sd_loc_utils.mod.F90
sd_loc_utils.mod.o:sd_loc_utils.mod.f90 cppt.mod fftmain_utils.mod \
                geq0mod.mod ions.mod kinds.mod parac.mod sfac.mod \
                system.mod

sd_nl2_utils.mod.f90:$(SRCDIR)/sd_nl2_utils.mod.F90
sd_nl2_utils.mod.o:sd_nl2_utils.mod.f90 elct.mod error_handling.mod \
                ions.mod kinds.mod nlps.mod parac.mod pslo.mod \
                sfac.mod sgpp.mod system.mod

sd_nl_utils.mod.f90:$(SRCDIR)/sd_nl_utils.mod.F90
sd_nl_utils.mod.o:sd_nl_utils.mod.f90 error_handling.mod ions.mod \
                kinds.mod nlps.mod parac.mod pslo.mod sgpp.mod \
                system.mod

sd_wannier_utils.mod.f90:$(SRCDIR)/sd_wannier_utils.mod.F90
sd_wannier_utils.mod.o:sd_wannier_utils.mod.f90 error_handling.mod \
                kinds.mod timer.mod utils.mod zeroing_utils.mod \
                kinds.mod error_handling.mod timer.mod wann.mod \
                kinds.mod error_handling.mod timer.mod wann.mod \
                zeroing_utils.mod kinds.mod error_handling.mod \
                timer.mod system.mod parac.mod wann.mod kinds.mod \
                error_handling.mod timer.mod prng_utils.mod \
                system.mod parac.mod wann.mod sd_wannier_utils.mod \
                utils.mod utils.mod

secder_driver.mod.f90:$(SRCDIR)/secder_driver.mod.F90
secder_driver.mod.o:secder_driver.mod.f90 adat.mod andp.mod \
                bsym.mod bsympnt.mod calc_alm_utils.mod cnst.mod \
                coor.mod copot_utils.mod cotr.mod ddip.mod \
                ddipo_utils.mod detdof_utils.mod dipo_utils.mod \
                dipomod.mod dynit_utils.mod ehpsi_utils.mod \
                elct.mod ener.mod error_handling.mod fileopen_utils.mod \
                fileopenmod.mod fint.mod forcedr_driver.mod \
                forcep_utils.mod gsize_utils.mod hessin_utils.mod \
                hessout_utils.mod initrun_driver.mod ions.mod \
                k_updwf_utils.mod kddipo_utils.mod kinds.mod \
                kpts.mod linres.mod lr_tddft_utils.mod lsforce_utils.mod \
                machine.mod mp_interface.mod nlcc.mod norm.mod \
                ortho_utils.mod parac.mod pbc_utils.mod phfac_utils.mod \
                poin.mod prop.mod pslo.mod rhoofr_c_utils.mod \
                rhoofr_utils.mod rmas.mod rnlsm_utils.mod ropt.mod \
                rrane_utils.mod secder_utils.mod setbsstate_utils.mod \
                setirec_utils.mod soft.mod spin.mod store_types.mod \
                symm.mod symtrz_utils.mod system.mod testex_utils.mod \
                updrho_utils.mod updwf_utils.mod utils.mod \
                wrener_utils.mod wrgeo_utils.mod wv30_utils.mod \
                zeroing_utils.mod

secder_utils.mod.f90:$(SRCDIR)/secder_utils.mod.F90
secder_utils.mod.o:secder_utils.mod.f90 adat.mod calc_alm_utils.mod \
                cell.mod cnst.mod copot_utils.mod ddipo_utils.mod \
                elct.mod error_handling.mod fileopen_utils.mod \
                fileopenmod.mod fint.mod forcedr_utils.mod \
                initrun_utils.mod ions.mod isos.mod kddipo_utils.mod \
                kinds.mod kpts.mod lr_tddft_utils.mod nlcc.mod \
                ortho_utils.mod parac.mod prop.mod pslo.mod \
                rhoofr_utils.mod rmas.mod rnlsm_utils.mod store_types.mod \
                symm.mod symtrz_utils.mod system.mod updrho_utils.mod \
                updwf_utils.mod utils.mod zeroing_utils.mod

secdpt_utils.mod.f90:$(SRCDIR)/secdpt_utils.mod.F90
secdpt_utils.mod.o:secdpt_utils.mod.f90 atwf.mod bsym.mod elct.mod \
                error_handling.mod fnlalloc_utils.mod ions.mod \
                kinds.mod linres.mod lr_in_utils.mod parac.mod \
                prmem_utils.mod prop.mod sdlinres_utils.mod \
                secder_driver.mod system.mod utils.mod vibana_utils.mod \
                zeroing_utils.mod

setbasis_utils.mod.f90:$(SRCDIR)/setbasis_utils.mod.F90
setbasis_utils.mod.o:setbasis_utils.mod.f90 adat.mod atom.mod \
                atwf.mod cdftmod.mod cnst.mod cppt.mod error_handling.mod \
                fitpack_utils.mod gvec.mod inscan_utils.mod \
                ions.mod kinds.mod lsfbtr_utils.mod mp_interface.mod \
                mw.mod parac.mod pimd.mod pslo.mod qspl.mod \
                readsr_utils.mod recpnew_utils.mod response_pmod.mod \
                sfac.mod sphe.mod system.mod timer.mod utils.mod \
                vdbp.mod ylmr_utils.mod zeroing_utils.mod

setbsstate_utils.mod.f90:$(SRCDIR)/setbsstate_utils.mod.F90
setbsstate_utils.mod.o:setbsstate_utils.mod.f90 bsym.mod elct.mod \
                hubbardu.mod kinds.mod mp_interface.mod parac.mod \
                ropt.mod spin.mod system.mod zeroing_utils.mod

setcnst_utils.mod.f90:$(SRCDIR)/setcnst_utils.mod.F90
setcnst_utils.mod.o:setcnst_utils.mod.f90 atoms_utils.mod cnst.mod \
                parac.mod soft.mod

set_cp_grp_utils.mod.f90:$(SRCDIR)/set_cp_grp_utils.mod.F90
set_cp_grp_utils.mod.o:set_cp_grp_utils.mod.f90 cp_grp_utils.mod \
                error_handling.mod kinds.mod machine.mod mp_interface.mod \
                parac.mod pimd.mod

setirec_utils.mod.f90:$(SRCDIR)/setirec_utils.mod.F90
setirec_utils.mod.o:setirec_utils.mod.f90 clas.mod cotr.mod \
                glemod.mod kinds.mod kpts.mod mp_interface.mod \
                nose.mod parac.mod shop_rest.mod store_types.mod \
                system.mod zeroing_utils.mod

setsc_utils.mod.f90:$(SRCDIR)/setsc_utils.mod.F90
setsc_utils.mod.o:setsc_utils.mod.f90 bc.mod cell.mod clas.mod \
                cnst.mod error_handling.mod gvec.mod isos.mod \
                kinds.mod latgen_utils.mod metr.mod molsym_utils.mod \
                mp_interface.mod parac.mod prcp.mod sort_utils.mod \
                symm.mod symmetry_utils.mod system.mod utils.mod \
                zeroing_utils.mod

setsys_utils.mod.f90:$(SRCDIR)/setsys_utils.mod.F90
setsys_utils.mod.o:setsys_utils.mod.f90 adat.mod atom.mod bsym.mod \
                cell.mod chksym_utils.mod clas.mod cnst.mod \
                cnst_dyn.mod coor.mod cotr.mod cplngs_utils.mod \
                dpot.mod elct.mod elct2.mod error_handling.mod \
                filnmod.mod fint.mod hubbardu.mod ions.mod \
                isos.mod kdpc.mod kinds.mod kpnt.mod kpts.mod \
                metr.mod mm_dim_utils.mod mm_dimmod.mod mm_extrap.mod \
                mm_input.mod mp_interface.mod multtb_utils.mod \
                nlcc.mod nlps.mod parac.mod pbc_utils.mod pslo.mod \
                ragg.mod response_pmod.mod rkpnt_utils.mod \
                rmas.mod ropt.mod rv30_utils.mod sfac.mod sgpp.mod \
                shock.mod shop.mod sphe.mod spin.mod store_types.mod \
                symm.mod system.mod vdbp.mod vdbt.mod velupi_utils.mod \
                wann.mod zeroing_utils.mod

sfac.mod.f90:   $(SRCDIR)/sfac.mod.F90
sfac.mod.o:     sfac.mod.f90 kinds.mod

sgpp.mod.f90:   $(SRCDIR)/sgpp.mod.F90
sgpp.mod.o:     sgpp.mod.f90 kinds.mod system.mod

shake_utils.mod.f90:$(SRCDIR)/shake_utils.mod.F90
shake_utils.mod.o:shake_utils.mod.f90 cnstfc_utils.mod cotr.mod \
                error_handling.mod ions.mod kinds.mod mm_dim_utils.mod \
                mm_dimmod.mod odiis_utils.mod parac.mod puttau_utils.mod \
                system.mod tpar.mod zeroing_utils.mod

shock.mod.f90:  $(SRCDIR)/shock.mod.F90
shock.mod.o:    shock.mod.f90 kinds.mod

shop_adds_utils.mod.f90:$(SRCDIR)/shop_adds_utils.mod.F90
shop_adds_utils.mod.o:shop_adds_utils.mod.f90 coor.mod elct.mod \
                ener.mod error_handling.mod fileopen_utils.mod \
                fileopenmod.mod geq0mod.mod kinds.mod machine.mod \
                mm_dimmod.mod mp_interface.mod parac.mod prng_utils.mod \
                shop.mod shop_ekinqm.mod shop_rest.mod spin.mod \
                system.mod zeroing_utils.mod

shop_ekinqm.mod.f90:$(SRCDIR)/shop_ekinqm.mod.F90
shop_ekinqm.mod.o:shop_ekinqm.mod.f90 kinds.mod

shop.mod.f90:   $(SRCDIR)/shop.mod.F90
shop.mod.o:     shop.mod.f90 kinds.mod

shop_rest_2.mod.f90:$(SRCDIR)/shop_rest_2.mod.F90
shop_rest_2.mod.o:shop_rest_2.mod.f90 kinds.mod

shop_rest.mod.f90:$(SRCDIR)/shop_rest.mod.F90
shop_rest.mod.o:shop_rest.mod.f90 kinds.mod

sh_tddft_utils.mod.f90:$(SRCDIR)/sh_tddft_utils.mod.F90
sh_tddft_utils.mod.o:sh_tddft_utils.mod.f90 adjmu_utils.mod \
                cnst.mod conv.mod coor.mod davidson_utils.mod \
                ddip.mod ddipo_utils.mod dotp_utils.mod elct.mod \
                ener.mod error_handling.mod fileopen_utils.mod \
                fileopenmod.mod fint.mod friesner_utils.mod \
                func.mod geq0mod.mod gsortho_utils.mod ions.mod \
                kinds.mod ksdiag_utils.mod linres.mod machine.mod \
                meta_multiple_walkers_utils.mod mm_dimmod.mod \
                mm_input.mod mp_interface.mod mw.mod opeigr_p_utils.mod \
                opeigr_utils.mod ovlap_utils.mod parac.mod \
                pimd.mod poin.mod prmem_utils.mod prng_utils.mod \
                randtowf_utils.mod readsr_utils.mod rhoofr_utils.mod \
                rmas.mod spin.mod symm.mod system.mod tauf.mod \
                td_nacvs_utils.mod timer.mod tpar.mod utils.mod \
                wrener_utils.mod zeroing_utils.mod

sh_utils.mod.f90:$(SRCDIR)/sh_utils.mod.F90
sh_utils.mod.o: sh_utils.mod.f90 cnst.mod error_handling.mod \
                kinds.mod linres.mod machine.mod mm_dimmod.mod \
                mm_input.mod mp_interface.mod nose.mod parac.mod \
                prng_utils.mod soc_types.mod system.mod zeroing_utils.mod

simple_model_p_utils.mod.f90:$(SRCDIR)/simple_model_p_utils.mod.F90
simple_model_p_utils.mod.o:simple_model_p_utils.mod.f90 error_handling.mod \
                kinds.mod mp_interface.mod ovlap_utils.mod \
                parac.mod system.mod timer.mod zeroing_utils.mod

simul.mod.f90:  $(SRCDIR)/simul.mod.F90
simul.mod.o:    simul.mod.f90 kinds.mod

sizeof_kinds.mod.f90:$(SRCDIR)/sizeof_kinds.mod.F90
sizeof_kinds.mod.o:sizeof_kinds.mod.f90 kinds.mod

soc.mod.f90:    $(SRCDIR)/soc.mod.F90
soc.mod.o:      soc.mod.f90 cnst.mod coor.mod cppt.mod dotp_utils.mod \
                elct.mod error_handling.mod fft_maxfft.mod \
                fftmain_utils.mod fftutil_utils.mod geq0mod.mod \
                ions.mod kinds.mod linres.mod lr_in_utils.mod \
                lr_tddft_utils.mod machine.mod mp_interface.mod \
                parac.mod prng_utils.mod response_pmod.mod \
                soc_types.mod specpt_utils.mod system.mod timer.mod \
                zeroing_utils.mod

soc_types.mod.f90:$(SRCDIR)/soc_types.mod.F90
soc_types.mod.o:soc_types.mod.f90 kinds.mod

softex_utils.mod.f90:$(SRCDIR)/softex_utils.mod.F90
softex_utils.mod.o:softex_utils.mod.f90 machine.mod parac.mod \
                soft.mod

soft.mod.f90:   $(SRCDIR)/soft.mod.F90
soft.mod.o:     soft.mod.f90

sort_utils.mod.f90:$(SRCDIR)/sort_utils.mod.F90
sort_utils.mod.o:sort_utils.mod.f90 error_handling.mod kinds.mod

special_functions.mod.f90:$(SRCDIR)/special_functions.mod.F90
special_functions.mod.o:special_functions.mod.f90 kinds.mod \
                kinds.mod kinds.mod kinds.mod kinds.mod kinds.mod \
                error_handling.mod kinds.mod kinds.mod error_handling.mod

specpt_utils.mod.f90:$(SRCDIR)/specpt_utils.mod.F90
specpt_utils.mod.o:specpt_utils.mod.f90 adjmu_utils.mod canon_utils.mod \
                conv.mod coor.mod corec_utils.mod cppt.mod \
                davidson_utils.mod ddip.mod ddipo_utils.mod \
                dftin_utils.mod dist_friesner_utils.mod dynit_utils.mod \
                efld.mod elct.mod ener.mod error_handling.mod \
                fftmain_utils.mod fint.mod fnlalloc_utils.mod \
                forcedr_driver.mod forcedr_utils.mod forcep_utils.mod \
                friesner_utils.mod func.mod gettrans_utils.mod \
                gsortho_utils.mod initrun_driver.mod initrun_utils.mod \
                ions.mod isos.mod kinds.mod kpts.mod ksdiag_utils.mod \
                linres.mod localize_utils.mod lr_diag_utils.mod \
                lr_in_utils.mod lr_tddft_drhoe.mod lr_tddft_utils.mod \
                lr_xcpot_utils.mod machine.mod mm_dim_utils.mod \
                mm_dimmod.mod mm_input.mod mm_qmmm_forcedr_utils.mod \
                mols.mod molstates_utils.mod mp_interface.mod \
                nlcc.mod norm.mod ortho_utils.mod ovlap_utils.mod \
                parac.mod pbc_utils.mod poin.mod prmem_utils.mod \
                pw_hfx_resp.mod pw_hfx_resp_types.mod readsr_utils.mod \
                reshaper.mod rho1pri_utils.mod rhoofr_utils.mod \
                rhopri_utils.mod rnlsm_utils.mod ropt.mod rv30_utils.mod \
                setirec_utils.mod soc_types.mod soft.mod sort_utils.mod \
                spin.mod stcop_utils.mod store_types.mod system.mod \
                tauf.mod td_force_utils.mod td_nacvs_utils.mod \
                td_os_utils.mod td_prop_utils.mod testex_utils.mod \
                timer.mod tpot.mod updwf_utils.mod utils.mod \
                wrener_utils.mod wv30_utils.mod zeroing_utils.mod

sphe.mod.f90:   $(SRCDIR)/sphe.mod.F90
sphe.mod.o:     sphe.mod.f90 kinds.mod

spin.mod.f90:   $(SRCDIR)/spin.mod.F90
spin.mod.o:     spin.mod.f90 kinds.mod

spsi_utils.mod.f90:$(SRCDIR)/spsi_utils.mod.F90
spsi_utils.mod.o:spsi_utils.mod.f90 cppt.mod cvan.mod error_handling.mod \
                ions.mod kinds.mod nlps.mod pslo.mod sfac.mod \
                system.mod timer.mod

ssic_utils.mod.f90:$(SRCDIR)/ssic_utils.mod.F90
ssic_utils.mod.o:ssic_utils.mod.f90 cppt.mod error_handling.mod \
                fft_maxfft.mod fftmain_utils.mod fftnew_utils.mod \
                geq0mod.mod kinds.mod mp_interface.mod parac.mod \
                system.mod zeroing_utils.mod

stagetrans_utils.mod.f90:$(SRCDIR)/stagetrans_utils.mod.F90
stagetrans_utils.mod.o:stagetrans_utils.mod.f90 ions.mod kinds.mod \
                pimd.mod

startpa_utils.mod.f90:$(SRCDIR)/startpa_utils.mod.F90
startpa_utils.mod.o:startpa_utils.mod.f90 error_handling.mod \
                fft_utils.mod geq0mod.mod machine.mod mp_interface.mod \
                parac.mod pimd_utils.mod readsr_utils.mod set_cp_grp_utils.mod \
                system.mod system_utils.mod utils.mod vdwcmod_utils.mod \
                zeroing_utils.mod

state_utils.mod.f90:$(SRCDIR)/state_utils.mod.F90
state_utils.mod.o:state_utils.mod.f90 cnst.mod cppt.mod fft.mod \
                geq0mod.mod kinds.mod system.mod

stcop_utils.mod.f90:$(SRCDIR)/stcop_utils.mod.F90
stcop_utils.mod.o:stcop_utils.mod.f90 dotp_utils.mod elct.mod \
                error_handling.mod hfx_drivers.mod kinds.mod \
                linres.mod mp_interface.mod parac.mod prng_utils.mod \
                rho1ofr_utils.mod spin.mod system.mod timer.mod \
                utils.mod v1ofrho1_utils.mod vpsi_utils.mod \
                zeroing_utils.mod

store_types.mod.f90:$(SRCDIR)/store_types.mod.F90
store_types.mod.o:store_types.mod.f90

str2.mod.f90:   $(SRCDIR)/str2.mod.F90
str2.mod.o:     str2.mod.f90 kinds.mod

stress_utils.mod.f90:$(SRCDIR)/stress_utils.mod.F90
stress_utils.mod.o:stress_utils.mod.f90 cnst.mod cppt.mod cvan.mod \
                elct.mod error_handling.mod ions.mod kdp.mod \
                kdp_stress_kin_utils.mod kdpc.mod kinds.mod \
                kpnt.mod kpts.mod metr.mod nlcc.mod nlps.mod \
                nlsm1_s_utils.mod parac.mod pbc_utils.mod prcp.mod \
                prmem_utils.mod pslo.mod putbet_utils.mod ragg.mod \
                rnlsm_utils.mod ropt.mod sfac.mod sgpp.mod \
                special_functions.mod spin.mod str2.mod strs.mod \
                system.mod timer.mod utils.mod vdwcmod.mod \
                zeroing_utils.mod

string_utils.mod.f90:$(SRCDIR)/string_utils.mod.F90
string_utils.mod.o:string_utils.mod.f90 kinds.mod

strs.mod.f90:   $(SRCDIR)/strs.mod.F90
strs.mod.o:     strs.mod.f90 kinds.mod

struc.mod.f90:  $(SRCDIR)/struc.mod.F90
struc.mod.o:    struc.mod.f90

struc_utils.mod.f90:$(SRCDIR)/struc_utils.mod.F90
struc_utils.mod.o:struc_utils.mod.f90 adat.mod cnst.mod coninp_utils.mod \
                constr_utils.mod empf.mod empfor_utils.mod \
                error_handling.mod fstart_utils.mod ions.mod \
                kinds.mod parac.mod struc.mod system.mod zeroing_utils.mod

sumfnl_utils.mod.f90:$(SRCDIR)/sumfnl_utils.mod.F90
sumfnl_utils.mod.o:sumfnl_utils.mod.f90 error_handling.mod \
                ions.mod kinds.mod mp_interface.mod nlps.mod \
                parac.mod sfac.mod system.mod timer.mod zeroing_utils.mod

summat_utils.mod.f90:$(SRCDIR)/summat_utils.mod.F90
summat_utils.mod.o:summat_utils.mod.f90 error_handling.mod \
                kinds.mod mp_interface.mod nlps.mod parac.mod \
                timer.mod kinds.mod error_handling.mod timer.mod \
                mp_interface.mod parac.mod

swap.mod.f90:   $(SRCDIR)/swap.mod.F90
swap.mod.o:     swap.mod.f90

symm4.mod.f90:  $(SRCDIR)/symm4.mod.F90
symm4.mod.o:    symm4.mod.f90 kinds.mod

symmetry_utils.mod.f90:$(SRCDIR)/symmetry_utils.mod.F90
symmetry_utils.mod.o:symmetry_utils.mod.f90 error_handling.mod \
                parac.mod

symm.mod.f90:   $(SRCDIR)/symm.mod.F90
symm.mod.o:     symm.mod.f90 kinds.mod system.mod

symtrz_utils.mod.f90:$(SRCDIR)/symtrz_utils.mod.F90
symtrz_utils.mod.o:symtrz_utils.mod.f90 cnst.mod cppt.mod error_handling.mod \
                fft_maxfft.mod fftmain_utils.mod fftnew_utils.mod \
                geq0mod.mod gvec.mod ions.mod kinds.mod metr.mod \
                mp_interface.mod parac.mod symm.mod symm4.mod \
                system.mod timer.mod zeroing_utils.mod

syscomb_utils.mod.f90:$(SRCDIR)/syscomb_utils.mod.F90
syscomb_utils.mod.o:syscomb_utils.mod.f90 cdftmod.mod cell.mod \
                cnst.mod coor.mod cppt.mod elct.mod error_handling.mod \
                fft.mod fft_maxfft.mod fftmain_utils.mod fftnew_utils.mod \
                filnmod.mod forcedr_driver.mod forcedr_utils.mod \
                geofile_utils.mod geq0mod.mod hpsi_utils.mod \
                ions.mod kinds.mod mp_interface.mod ortho_utils.mod \
                ovlap_utils.mod parac.mod prden.mod rhoofr_utils.mod \
                rhopri_utils.mod rmas.mod rv30_utils.mod setirec_utils.mod \
                spin.mod store_types.mod system.mod utils.mod \
                wv30_utils.mod zeroing_utils.mod

sysdepend.o:    $(SRCDIR)/sysdepend.c

sysin_utils.mod.f90:$(SRCDIR)/sysin_utils.mod.F90
sysin_utils.mod.o:sysin_utils.mod.f90 cdftmod.mod cell.mod \
                clas.mod cnst.mod cplngs_utils.mod cplngsmod.mod \
                ddip.mod dg.mod elct.mod elct2.mod error_handling.mod \
                fcas.mod gvec.mod hfxmod.mod inscan_utils.mod \
                isos.mod kdpc.mod kinds.mod kpnt.mod kpts.mod \
                mp_interface.mod parac.mod prcp.mod readsr_utils.mod \
                ropt.mod shock.mod sphe.mod spin.mod symm.mod \
                system.mod zeroing_utils.mod

system.mod.f90: $(SRCDIR)/system.mod.F90
system.mod.o:   system.mod.f90 kinds.mod

system_utils.mod.f90:$(SRCDIR)/system_utils.mod.F90
system_utils.mod.o:system_utils.mod.f90 error_handling.mod \
                system.mod

tauf.mod.f90:   $(SRCDIR)/tauf.mod.F90
tauf.mod.o:     tauf.mod.f90 kinds.mod

tauofr_utils.mod.f90:$(SRCDIR)/tauofr_utils.mod.F90
tauofr_utils.mod.o:tauofr_utils.mod.f90 cnst.mod cppt.mod elct.mod \
                error_handling.mod fft_maxfft.mod fftmain_utils.mod \
                geq0mod.mod kinds.mod parac.mod pslo.mod spin.mod \
                system.mod tauf.mod timer.mod zeroing_utils.mod

tbxc.mod.f90:   $(SRCDIR)/tbxc.mod.F90
tbxc.mod.o:     tbxc.mod.f90 kinds.mod

td_cayley_utils.mod.f90:$(SRCDIR)/td_cayley_utils.mod.F90
td_cayley_utils.mod.o:td_cayley_utils.mod.f90 cnst.mod coor.mod \
                dipo_utils.mod dipomod.mod eicalc_utils.mod \
                ener.mod error_handling.mod fft_maxfft.mod \
                fftmain_utils.mod geq0mod.mod hpsi_utils.mod \
                ions.mod kinds.mod mp_interface.mod parac.mod \
                rhoofr_c_utils.mod soft.mod spin.mod system.mod \
                td_input.mod td_utils.mod testex_utils.mod \
                timer.mod utils.mod vofrho_utils.mod zeroing_utils.mod

td_dav_utils.mod.f90:$(SRCDIR)/td_dav_utils.mod.F90
td_dav_utils.mod.o:td_dav_utils.mod.f90 cnst.mod dotp_utils.mod \
                elct.mod error_handling.mod kinds.mod ksdiag_utils.mod \
                linres.mod lr_force_utils.mod lr_ortho_utils.mod \
                machine.mod mp_interface.mod parac.mod pw_hfx_resp_types.mod \
                randtowf_utils.mod soft.mod sort_utils.mod \
                system.mod testex_utils.mod timer.mod zeroing_utils.mod

td_force_utils.mod.f90:$(SRCDIR)/td_force_utils.mod.F90
td_force_utils.mod.o:td_force_utils.mod.f90 afbdr_utils.mod \
                cppt.mod elct.mod error_handling.mod fft.mod \
                fftmain_utils.mod ions.mod kinds.mod ksdiag_utils.mod \
                linres.mod lr_ortho_utils.mod lr_xcpot_utils.mod \
                machine.mod mp_interface.mod opt_lr_utils.mod \
                ovlap_utils.mod parac.mod poin.mod rho1ofr_utils.mod \
                rnlfor_utils.mod rnlsm_utils.mod ropt.mod rotate_utils.mod \
                sfac.mod spin.mod system.mod timer.mod v1ofrho1_utils.mod \
                vpsi_utils.mod vtd2_utils.mod wrgeo_utils.mod \
                zeroing_utils.mod

td_input.mod.f90:$(SRCDIR)/td_input.mod.F90
td_input.mod.o: td_input.mod.f90 kinds.mod

td_mm_qmmm_forcedr_utils.mod.f90:$(SRCDIR)/td_mm_qmmm_forcedr_utils.mod.F90
td_mm_qmmm_forcedr_utils.mod.o:td_mm_qmmm_forcedr_utils.mod.f90 \
                ener.mod error_handling.mod forcedr_driver.mod \
                kinds.mod machine.mod mm_dim_utils.mod mm_dimmod.mod \
                mm_input.mod mp_interface.mod parac.mod puttau_utils.mod \
                rhoofr_utils.mod system.mod zeroing_utils.mod

td_nacvs_utils.mod.f90:$(SRCDIR)/td_nacvs_utils.mod.F90
td_nacvs_utils.mod.o:td_nacvs_utils.mod.f90 canon_utils.mod \
                cell.mod cnst.mod constr_utils.mod coor.mod \
                cplngs_utils.mod cplngsmod.mod cppt.mod dotp_utils.mod \
                eicalc_utils.mod eind_ii_utils.mod eind_loc_utils.mod \
                eind_nl_utils.mod ekinpp_utils.mod elct.mod \
                error_handling.mod fftmain_utils.mod fileopen_utils.mod \
                fileopenmod.mod fnlalloc_utils.mod forcedr_driver.mod \
                forcedr_utils.mod forcep_utils.mod geq0mod.mod \
                ions.mod kinds.mod ksdiag_utils.mod linres.mod \
                lr_tddft_drhoe.mod lr_xcpot_utils.mod mp_interface.mod \
                nl_res_utils.mod nlps.mod nose.mod opt_lr_utils.mod \
                parac.mod poin.mod pslo.mod rho1ofr_utils.mod \
                rhoofr_utils.mod rmas.mod rnlsm1_utils.mod \
                rnlsm2_utils.mod rnlsm_2d_utils.mod rnlsm_utils.mod \
                ropt.mod sfac.mod sgpp.mod spin.mod symm.mod \
                symtrz_utils.mod system.mod td_force_utils.mod \
                timer.mod v1ofrho1_utils.mod zeroing_utils.mod

td_nhdav_utils.mod.f90:$(SRCDIR)/td_nhdav_utils.mod.F90
td_nhdav_utils.mod.o:td_nhdav_utils.mod.f90 elct.mod error_handling.mod \
                kinds.mod ksdiag_utils.mod linres.mod lr_force_utils.mod \
                lr_ortho_utils.mod machine.mod mp_interface.mod \
                parac.mod soft.mod sort_utils.mod system.mod \
                td_dav_utils.mod testex_utils.mod timer.mod \
                zeroing_utils.mod

tdnlfor_utils.mod.f90:$(SRCDIR)/tdnlfor_utils.mod.F90
tdnlfor_utils.mod.o:tdnlfor_utils.mod.f90 error_handling.mod \
                ions.mod kinds.mod nlps.mod parac.mod pslo.mod \
                sfac.mod sgpp.mod system.mod timer.mod

td_os_berry_utils.mod.f90:$(SRCDIR)/td_os_berry_utils.mod.F90
td_os_berry_utils.mod.o:td_os_berry_utils.mod.f90 cnst.mod \
                ddip.mod ddipo_utils.mod error_handling.mod \
                geq0mod.mod kinds.mod mp_interface.mod opeigr_p_utils.mod \
                parac.mod system.mod timer.mod utils.mod zeroing_utils.mod

td_os_utils.mod.f90:$(SRCDIR)/td_os_utils.mod.F90
td_os_utils.mod.o:td_os_utils.mod.f90 error_handling.mod kinds.mod \
                linres.mod opeigr_utils.mod parac.mod prmem_utils.mod \
                symm.mod system.mod td_os_berry_utils.mod

td_pcg_utils.mod.f90:$(SRCDIR)/td_pcg_utils.mod.F90
td_pcg_utils.mod.o:td_pcg_utils.mod.f90 cnst.mod elct.mod error_handling.mod \
                kinds.mod ksdiag_utils.mod linres.mod lr_force_utils.mod \
                lr_ortho_utils.mod machine.mod mp_interface.mod \
                orbrot_utils.mod parac.mod soft.mod system.mod \
                td_dav_utils.mod timer.mod zeroing_utils.mod

td_prop_utils.mod.f90:$(SRCDIR)/td_prop_utils.mod.F90
td_prop_utils.mod.o:td_prop_utils.mod.f90 atomc_utils.mod cppt.mod \
                eicalc_utils.mod elstpo_utils.mod error_handling.mod \
                espchg_utils.mod fft_maxfft.mod fftmain_utils.mod \
                fftnew_utils.mod geq0mod.mod ions.mod isos.mod \
                kinds.mod parac.mod system.mod timer.mod zeroing_utils.mod \
                kinds.mod error_handling.mod timer.mod mp_interface.mod \
                fft_maxfft.mod system.mod parac.mod cppt.mod \
                ions.mod isos.mod cnst.mod geq0mod.mod fftmain_utils.mod

td_utils.mod.f90:$(SRCDIR)/td_utils.mod.F90
td_utils.mod.o: td_utils.mod.f90 calc_pij_utils.mod cnst.mod \
                coor.mod cppt.mod densrd_utils.mod dg.mod dipo_utils.mod \
                dipomod.mod efld.mod eicalc_utils.mod elct.mod \
                ener.mod error_handling.mod fft.mod fft_maxfft.mod \
                fftmain_utils.mod fftnew_utils.mod fileopen_utils.mod \
                fileopenmod.mod geq0mod.mod gvec.mod hpsi_utils.mod \
                ions.mod kinds.mod kpclean_utils.mod kpnt.mod \
                machine.mod mp_interface.mod parac.mod phfac_utils.mod \
                readsr_utils.mod rhoofr_c_utils.mod rhopri_utils.mod \
                rnlrh_utils.mod rnlsm_utils.mod ropt.mod spin.mod \
                system.mod td_input.mod utils.mod vofrho_utils.mod \
                zeroing_utils.mod

temps.mod.f90:  $(SRCDIR)/temps.mod.F90
temps.mod.o:    temps.mod.f90 kinds.mod

testex_utils.mod.f90:$(SRCDIR)/testex_utils.mod.F90
testex_utils.mod.o:testex_utils.mod.f90 error_handling.mod \
                fileopenmod.mod kinds.mod mp_interface.mod \
                mw.mod parac.mod pimd.mod softex_utils.mod \
                system.mod timer.mod

teststore_utils.mod.f90:$(SRCDIR)/teststore_utils.mod.F90
teststore_utils.mod.o:teststore_utils.mod.f90 store_types.mod \
                system.mod

thread_view_types.mod.f90:$(SRCDIR)/thread_view_types.mod.F90
thread_view_types.mod.o:thread_view_types.mod.f90 error_handling.mod

thread_view_utils.mod.f90:$(SRCDIR)/thread_view_utils.mod.F90
thread_view_utils.mod.o:thread_view_utils.mod.f90 error_handling.mod \
                thread_view_types.mod

time.mod.f90:   $(SRCDIR)/time.mod.F90
time.mod.o:     time.mod.f90 kinds.mod

timer.mod.f90:  $(SRCDIR)/timer.mod.F90
timer.mod.o:    timer.mod.f90 envj.mod kinds.mod machine.mod \
                mp_interface.mod parac.mod system.mod time.mod \
                zeroing_utils.mod kinds.mod machine.mod time.mod

timetag.f90:    $(SRCDIR)/timetag.F90
timetag.o:      timetag.f90 kinds.mod error_handling.mod parac.mod

totstr_utils.mod.f90:$(SRCDIR)/totstr_utils.mod.F90
totstr_utils.mod.o:totstr_utils.mod.f90 cnst.mod cnst_dyn.mod \
                kinds.mod metr.mod mp_interface.mod parac.mod \
                prcp.mod strs.mod symtrz_utils.mod system.mod \
                vdwcmod.mod

tpar.mod.f90:   $(SRCDIR)/tpar.mod.F90
tpar.mod.o:     tpar.mod.f90 kinds.mod system.mod

tpot.mod.f90:   $(SRCDIR)/tpot.mod.F90
tpot.mod.o:     tpot.mod.f90 kinds.mod

transme.mod.f90:$(SRCDIR)/transme.mod.F90
transme.mod.o:  transme.mod.f90 kinds.mod

transme_utils.mod.f90:$(SRCDIR)/transme_utils.mod.F90
transme_utils.mod.o:transme_utils.mod.f90 cdftmod.mod cp_grp_utils.mod \
                cppt.mod elct.mod error_handling.mod filnmod.mod \
                geq0mod.mod gvec.mod ions.mod kinds.mod mp_interface.mod \
                parac.mod part_1d.mod setirec_utils.mod spin.mod \
                store_types.mod system.mod timer.mod transmemod.mod \
                zeroing_utils.mod

tst2min_inp_utils.mod.f90:$(SRCDIR)/tst2min_inp_utils.mod.F90
tst2min_inp_utils.mod.o:tst2min_inp_utils.mod.f90 cnst_dyn.mod \
                error_handling.mod ions.mod kinds.mod mm_dim_utils.mod \
                mm_dimmod.mod mm_input.mod parac.mod readsr_utils.mod \
                zeroing_utils.mod

tst2min_utils.mod.f90:$(SRCDIR)/tst2min_utils.mod.F90
tst2min_utils.mod.o:tst2min_utils.mod.f90 cnst_dyn.mod fileopen_utils.mod \
                fileopenmod.mod kinds.mod meta_colvar_inp_utils.mod \
                mp_interface.mod parac.mod readsr_utils.mod \
                soft.mod

up3_p_utils.mod.f90:$(SRCDIR)/up3_p_utils.mod.F90
up3_p_utils.mod.o:up3_p_utils.mod.f90 cnst.mod elct.mod error_handling.mod \
                geq0mod.mod gsortho_utils.mod kinds.mod kpnt.mod \
                kpts.mod mp_interface.mod ortho_utils.mod parac.mod \
                rkpnt_utils.mod rnlsm_utils.mod setirec_utils.mod \
                sphe.mod system.mod utils.mod wv30_utils.mod \
                zeroing_utils.mod

updrho_utils.mod.f90:$(SRCDIR)/updrho_utils.mod.F90
updrho_utils.mod.o:updrho_utils.mod.f90 adjmu_utils.mod andp.mod \
                bogol_utils.mod broy.mod conv.mod davidson_utils.mod \
                dist_friesner_utils.mod ehpsi_utils.mod elct.mod \
                elct2.mod ener.mod error_handling.mod fint.mod \
                friesner_c_p_utils.mod friesner_c_utils.mod \
                friesner_utils.mod frsblk_c_utils.mod frsblk_utils.mod \
                func.mod hfx_drivers.mod hpsi_utils.mod hubbardu.mod \
                hubbardu_utils.mod k_diis_rhofix_utils.mod \
                k_forces_utils.mod kdp.mod kdp_diag_utils.mod \
                kdp_prep_utils.mod kdp_rho_utils.mod kdpc.mod \
                kinds.mod kpnt.mod kpts.mod ksdiag_utils.mod \
                mixing_g_utils.mod mixing_r_utils.mod mp_interface.mod \
                norm.mod parac.mod pcgrad_driver.mod pslo.mod \
                ptheory_utils.mod rhoofr_c_utils.mod rhoofr_utils.mod \
                rnlfor_utils.mod rnlrh_utils.mod rnlsm_utils.mod \
                ropt.mod rpiiint_utils.mod soft.mod spin.mod \
                stress_utils.mod symtrz_utils.mod system.mod \
                tauf.mod testex_utils.mod timer.mod utils.mod \
                vbeta_utils.mod vdw_utils.mod vdwcmod.mod vofrho_utils.mod \
                zeroing_utils.mod

updwf_p_utils.mod.f90:$(SRCDIR)/updwf_p_utils.mod.F90
updwf_p_utils.mod.o:updwf_p_utils.mod.f90 csize_utils.mod error_handling.mod \
                forces_p_utils.mod geq0mod.mod kinds.mod norm.mod \
                odiis_p_utils.mod ortho_utils.mod pcgrad_p_utils.mod \
                perturbation_p_utils.mod response_pmod.mod \
                ropt.mod simple_model_p_utils.mod soft.mod \
                system.mod testex_utils.mod tpar.mod utils.mod

updwf_utils.mod.f90:$(SRCDIR)/updwf_utils.mod.F90
updwf_utils.mod.o:updwf_utils.mod.f90 adapttol_utils.mod forcedr_driver.mod \
                forcedr_utils.mod geq0mod.mod hesele_utils.mod \
                hubbardu.mod kinds.mod mm_dim_utils.mod mm_dimmod.mod \
                mm_input.mod mm_qmmm_forcedr_utils.mod norm.mod \
                ortho_utils.mod parac.mod pcgrad_driver.mod \
                pcgrad_utils.mod pslo.mod rnlsm_utils.mod ropt.mod \
                soft.mod system.mod testex_utils.mod timer.mod \
                tpar.mod utils.mod zeroing_utils.mod

util_p_utils.mod.f90:$(SRCDIR)/util_p_utils.mod.F90
util_p_utils.mod.o:util_p_utils.mod.f90 cppt.mod error_handling.mod \
                fileopen_utils.mod kinds.mod parac.mod system.mod \
                mp_interface.mod readsr_utils.mod system.mod \
                parac.mod ions.mod coor.mod metr.mod response_pmod.mod \
                fileopen_utils.mod fileopen_utils.mod fileopenmod.mod \
                timer.mod system.mod parac.mod cppt.mod sfac.mod \
                geq0mod.mod timer.mod system.mod parac.mod \
                cppt.mod sfac.mod geq0mod.mod

utils.mod.f90:  $(SRCDIR)/utils.mod.F90
utils.mod.o:    utils.mod.f90 error_handling.mod geq0mod.mod \
                kinds.mod parac.mod reshaper.mod system.mod \
                zeroing_utils.mod kinds.mod kinds.mod kinds.mod \
                kinds.mod kinds.mod kinds.mod

u_upd_exp_sum_utils.mod.f90:$(SRCDIR)/u_upd_exp_sum_utils.mod.F90
u_upd_exp_sum_utils.mod.o:u_upd_exp_sum_utils.mod.f90 error_handling.mod \
                g_loc.mod g_loc_optim_utils.mod g_loc_util_utils.mod \
                kinds.mod parac.mod system.mod timer.mod utils.mod \
                zeroing_utils.mod kinds.mod error_handling.mod \
                timer.mod system.mod parac.mod

u_upd_exp_utils.mod.f90:$(SRCDIR)/u_upd_exp_utils.mod.F90
u_upd_exp_utils.mod.o:u_upd_exp_utils.mod.f90 error_handling.mod \
                g_loc.mod g_loc_util_utils.mod kinds.mod parac.mod \
                system.mod timer.mod utils.mod zeroing_utils.mod \
                znum_mat_utils.mod

u_upd_spread_sum_utils.mod.f90:$(SRCDIR)/u_upd_spread_sum_utils.mod.F90
u_upd_spread_sum_utils.mod.o:u_upd_spread_sum_utils.mod.f90 \
                error_handling.mod g_loc.mod g_loc_util_utils.mod \
                kinds.mod parac.mod system.mod timer.mod utils.mod \
                zeroing_utils.mod

v1ofrho1_utils.mod.f90:$(SRCDIR)/v1ofrho1_utils.mod.F90
v1ofrho1_utils.mod.o:v1ofrho1_utils.mod.f90 cppt.mod dd_xc_ana_utils.mod \
                dd_xc_utils.mod error_handling.mod fftmain_utils.mod \
                fftnew_utils.mod geq0mod.mod isos.mod kinds.mod \
                linres.mod parac.mod poin.mod spin.mod system.mod \
                timer.mod zeroing_utils.mod

v1ofrho_p_utils.mod.f90:$(SRCDIR)/v1ofrho_p_utils.mod.F90
v1ofrho_p_utils.mod.o:v1ofrho_p_utils.mod.f90 cppt.mod error_handling.mod \
                fft_maxfft.mod geq0mod.mod isos.mod kinds.mod \
                reshaper.mod response_pmod.mod spin.mod system.mod \
                timer.mod v1xc_p_utils.mod zeroing_utils.mod

v1xc_p_utils.mod.f90:$(SRCDIR)/v1xc_p_utils.mod.F90
v1xc_p_utils.mod.o:v1xc_p_utils.mod.f90 cppt.mod error_handling.mod \
                fft_maxfft.mod fftmain_utils.mod func.mod gcener_utils.mod \
                graden_utils.mod kinds.mod nlcc.mod parac.mod \
                reshaper.mod response_pmod.mod spin.mod system.mod \
                timer.mod utils.mod xcener_utils.mod zeroing_utils.mod

vbeta_utils.mod.f90:$(SRCDIR)/vbeta_utils.mod.F90
vbeta_utils.mod.o:vbeta_utils.mod.f90 cnst.mod cppt.mod fft_maxfft.mod \
                fftmain_utils.mod fint.mod geq0mod.mod kinds.mod \
                kpnt.mod kpts.mod nlps.mod parac.mod sfac.mod \
                system.mod timer.mod zeroing_utils.mod

vdbinit_utils.mod.f90:$(SRCDIR)/vdbinit_utils.mod.F90
vdbinit_utils.mod.o:vdbinit_utils.mod.f90 aavan.mod cnst.mod \
                cppt.mod cvan.mod error_handling.mod fitpack_utils.mod \
                ions.mod kinds.mod mp_interface.mod nlps.mod \
                parac.mod pslo.mod qspl.mod qvan1_utils.mod \
                radin_utils.mod system.mod timer.mod vdbp.mod \
                ylmr_utils.mod zeroing_utils.mod

vdbp.mod.f90:   $(SRCDIR)/vdbp.mod.F90
vdbp.mod.o:     vdbp.mod.f90 kinds.mod system.mod

vdbt.mod.f90:   $(SRCDIR)/vdbt.mod.F90
vdbt.mod.o:     vdbt.mod.f90 system.mod

vdwcmod.mod.f90:$(SRCDIR)/vdwcmod.mod.F90
vdwcmod.mod.o:  vdwcmod.mod.f90 kinds.mod

vdwcmod_utils.mod.f90:$(SRCDIR)/vdwcmod_utils.mod.F90
vdwcmod_utils.mod.o:vdwcmod_utils.mod.f90 error_handling.mod \
                vdwcmod.mod

vdwin_utils.mod.f90:$(SRCDIR)/vdwin_utils.mod.F90
vdwin_utils.mod.o:vdwin_utils.mod.f90 adat.mod cnst.mod dcacp_utils.mod \
                error_handling.mod inscan_utils.mod ions.mod \
                kinds.mod mp_interface.mod parac.mod readsr_utils.mod \
                system.mod timer.mod vdwcmod.mod wann.mod

vdw_utils.mod.f90:$(SRCDIR)/vdw_utils.mod.F90
vdw_utils.mod.o:vdw_utils.mod.f90 adat.mod cnst.mod error_handling.mod \
                fileopen_utils.mod fileopenmod.mod ions.mod \
                kinds.mod meta_multiple_walkers_utils.mod metr.mod \
                mp_interface.mod mw.mod parac.mod pbc_utils.mod \
                pimd.mod sort_utils.mod strs.mod system.mod \
                timer.mod vdwcmod.mod wrgeo_utils.mod zeroing_utils.mod

vdw_wf_alloc_utils.mod.f90:$(SRCDIR)/vdw_wf_alloc_utils.mod.F90
vdw_wf_alloc_utils.mod.o:vdw_wf_alloc_utils.mod.f90 elct.mod \
                error_handling.mod ions.mod pimd.mod system.mod \
                vdwcmod.mod

velocitinp_utils.mod.f90:$(SRCDIR)/velocitinp_utils.mod.F90
velocitinp_utils.mod.o:velocitinp_utils.mod.f90 coor.mod error_handling.mod \
                ions.mod kinds.mod parac.mod readsr_utils.mod \
                system.mod

velupa_utils.mod.f90:$(SRCDIR)/velupa_utils.mod.F90
velupa_utils.mod.o:velupa_utils.mod.f90 harm.mod kinds.mod \
                system.mod timer.mod tpar.mod

velupi_utils.mod.f90:$(SRCDIR)/velupi_utils.mod.F90
velupi_utils.mod.o:velupi_utils.mod.f90 cnst.mod ions.mod kinds.mod \
                metr.mod parac.mod pimd.mod prcp.mod puttau_utils.mod \
                rmas.mod system.mod tpar.mod utils.mod zeroing_utils.mod

vepsup_utils.mod.f90:$(SRCDIR)/vepsup_utils.mod.F90
vepsup_utils.mod.o:vepsup_utils.mod.f90 cnst.mod ions.mod kinds.mod \
                metr.mod rmas.mod shock.mod system.mod tpar.mod \
                zeroing_utils.mod

vgsortho_utils.mod.f90:$(SRCDIR)/vgsortho_utils.mod.F90
vgsortho_utils.mod.o:vgsortho_utils.mod.f90 dotp_utils.mod \
                geq0mod.mod kinds.mod mp_interface.mod parac.mod \
                timer.mod

vhk_utils.mod.f90:$(SRCDIR)/vhk_utils.mod.F90
vhk_utils.mod.o:vhk_utils.mod.f90 error_handling.mod fft_maxfft.mod \
                gcener_utils.mod graden_utils.mod kinds.mod \
                nlcc.mod spin.mod system.mod timer.mod xcener_utils.mod \
                zeroing_utils.mod

vibana_utils.mod.f90:$(SRCDIR)/vibana_utils.mod.F90
vibana_utils.mod.o:vibana_utils.mod.f90 adat.mod coor.mod cotr.mod \
                detdof_utils.mod error_handling.mod forcep_utils.mod \
                hessin_utils.mod hessout_utils.mod initrun_driver.mod \
                initrun_utils.mod ions.mod kinds.mod lscal.mod \
                parac.mod rmas.mod secder_utils.mod store_types.mod \
                symm.mod symtrz_utils.mod system.mod utils.mod \
                zeroing_utils.mod

vlocst_utils.mod.f90:$(SRCDIR)/vlocst_utils.mod.F90
vlocst_utils.mod.o:vlocst_utils.mod.f90 cppt.mod ffsum_utils.mod \
                fft_maxfft.mod geq0mod.mod kinds.mod pslo.mod \
                str2.mod strs.mod system.mod timer.mod zeroing_utils.mod

voa_p_utils.mod.f90:$(SRCDIR)/voa_p_utils.mod.F90
voa_p_utils.mod.o:voa_p_utils.mod.f90 adat.mod cnst.mod coor.mod \
                ddip.mod dipo_utils.mod dipomod.mod dotp_utils.mod \
                eicalc_utils.mod elct.mod error_handling.mod \
                fftmain_utils.mod fileopen_utils.mod fileopenmod.mod \
                fnonloc_p_utils.mod fnonloc_utils.mod h0psi1_p_utils.mod \
                ions.mod kinds.mod localize_utils.mod mp_interface.mod \
                nlps.mod nmr_position_p_utils.mod nmr_util_p_utils.mod \
                parac.mod prop.mod response_pmod.mod rhoofr_p_utils.mod \
                rhoofr_utils.mod rmas.mod rnlsm_utils.mod ropt.mod \
                rotate_utils.mod rwfopt_p_utils.mod sfac.mod \
                spin.mod system.mod timer.mod utils.mod wann.mod \
                zeroing_utils.mod

vofrhoa_utils.mod.f90:$(SRCDIR)/vofrhoa_utils.mod.F90
vofrhoa_utils.mod.o:vofrhoa_utils.mod.f90 eextern_utils.mod \
                efld.mod eicalc_utils.mod elct.mod ener.mod \
                error_handling.mod fftmain_utils.mod htrstr_utils.mod \
                kinds.mod mp_interface.mod nvtx_utils.mod parac.mod \
                potfor_utils.mod ppener_utils.mod simulmod.mod \
                state_utils.mod system.mod td_input.mod timer.mod \
                vlocst_utils.mod zeroing_utils.mod

vofrhob_utils.mod.f90:$(SRCDIR)/vofrhob_utils.mod.F90
vofrhob_utils.mod.o:vofrhob_utils.mod.f90 cofor_utils.mod corec_utils.mod \
                cppt.mod cvan.mod efld.mod elct.mod ener.mod \
                error_handling.mod extpotmod.mod fft_maxfft.mod \
                fftmain_utils.mod fftnew_utils.mod gcener_utils.mod \
                geq0mod.mod graden_utils.mod isos.mod kinds.mod \
                linres.mod meta_localizespin_utils.mod newd_utils.mod \
                nlcc.mod nlccstr_utils.mod nlps.mod nvtx_utils.mod \
                parac.mod pslo.mod reshaper.mod saop_utils.mod \
                sfac.mod spin.mod ssic_utils.mod str2.mod strs.mod \
                system.mod tauf.mod timer.mod tpot.mod utils.mod \
                xcener_utils.mod zeroing_utils.mod

vofrhoc_utils.mod.f90:$(SRCDIR)/vofrhoc_utils.mod.F90
vofrhoc_utils.mod.o:vofrhoc_utils.mod.f90 cnst.mod cofor_utils.mod \
                cppt.mod eextern_utils.mod efld.mod eicalc_utils.mod \
                elct.mod ener.mod error_handling.mod extpotmod.mod \
                fcas.mod fftmain_utils.mod fftnew_utils.mod \
                gcener_utils.mod geq0mod.mod graden_utils.mod \
                ions.mod kinds.mod mp_interface.mod nlcc.mod \
                parac.mod potfor_utils.mod ppener_utils.mod \
                rnlfor_utils.mod ropt.mod rpiiint_utils.mod \
                sfac.mod spin.mod system.mod timer.mod utils.mod \
                xcener_utils.mod zeroing_utils.mod

vofrhoh_utils.mod.f90:$(SRCDIR)/vofrhoh_utils.mod.F90
vofrhoh_utils.mod.o:vofrhoh_utils.mod.f90 cppt.mod eextern_utils.mod \
                efld.mod eicalc_utils.mod ener.mod error_handling.mod \
                fftmain_utils.mod fftnew_utils.mod geq0mod.mod \
                hip_utils.mod ions.mod isos.mod kinds.mod mp_interface.mod \
                parac.mod sfac.mod simulmod.mod system.mod \
                td_input.mod timer.mod

vofrhos_utils.mod.f90:$(SRCDIR)/vofrhos_utils.mod.F90
vofrhos_utils.mod.o:vofrhos_utils.mod.f90 cnst.mod cofor_utils.mod \
                corec_utils.mod cppt.mod efld.mod ener.mod \
                error_handling.mod extpotmod.mod fftmain_utils.mod \
                fftnew_utils.mod func.mod gcener_utils.mod \
                geq0mod.mod graden_utils.mod isos.mod kinds.mod \
                nlcc.mod parac.mod potfor_utils.mod pslo.mod \
                spin.mod str2.mod system.mod timer.mod utils.mod \
                vofrhoc_utils.mod xcener_utils.mod zeroing_utils.mod

vofrhot_utils.mod.f90:$(SRCDIR)/vofrhot_utils.mod.F90
vofrhot_utils.mod.o:vofrhot_utils.mod.f90 eextern_utils.mod \
                efld.mod eicalc_utils.mod elct.mod ener.mod \
                error_handling.mod fftmain_utils.mod kinds.mod \
                mp_interface.mod parac.mod potfor_utils.mod \
                ppener_utils.mod simulmod.mod system.mod td_input.mod \
                timer.mod

vofrho_utils.mod.f90:$(SRCDIR)/vofrho_utils.mod.F90
vofrho_utils.mod.o:vofrho_utils.mod.f90 error_handling.mod \
                forcep_utils.mod isos.mod kinds.mod poin.mod \
                spin.mod store_types.mod system.mod timer.mod \
                vofrhoa_utils.mod vofrhob_utils.mod vofrhoc_utils.mod \
                vofrhoh_utils.mod vofrhos_utils.mod vofrhot_utils.mod

vpsi_lse_utils.mod.f90:$(SRCDIR)/vpsi_lse_utils.mod.F90
vpsi_lse_utils.mod.o:vpsi_lse_utils.mod.f90 cnst.mod cppt.mod \
                dotp_utils.mod ener.mod error_handling.mod \
                fft_maxfft.mod fftmain_utils.mod geq0mod.mod \
                kinds.mod mp_interface.mod parac.mod spin.mod \
                system.mod timer.mod zeroing_utils.mod

vpsi_p_utils.mod.f90:$(SRCDIR)/vpsi_p_utils.mod.F90
vpsi_p_utils.mod.o:vpsi_p_utils.mod.f90 cnst.mod cppt.mod elct.mod \
                error_handling.mod fft_maxfft.mod fftmain_utils.mod \
                geq0mod.mod kinds.mod kpts.mod parac.mod prcp.mod \
                reshaper.mod response_pmod.mod spin.mod system.mod \
                timer.mod zeroing_utils.mod

vpsi_utils.mod.f90:$(SRCDIR)/vpsi_utils.mod.F90
vpsi_utils.mod.o:vpsi_utils.mod.f90 cdftmod.mod cnst.mod cp_cuda_types.mod \
                cp_cufft_types.mod cp_cuvpsi_types.mod cp_cuvpsi_utils.mod \
                cp_cuwfn_types.mod cp_grp_utils.mod cppt.mod \
                cuda_types.mod cuda_utils.mod cuuser_utils.mod \
                dg.mod efld.mod error_handling.mod fft.mod \
                fftmain_utils.mod fftnew_utils.mod geq0mod.mod \
                kinds.mod kpclean_utils.mod kpnt.mod kpts.mod \
                lr_xcpot_utils.mod mp_interface.mod parac.mod \
                part_1d.mod prcp.mod reshaper.mod rswfmod.mod \
                special_functions.mod spin.mod state_utils.mod \
                system.mod td_input.mod thread_view_types.mod \
                thread_view_utils.mod timer.mod vpsi_lse_utils.mod \
                vtaupsi_utils.mod zeroing_utils.mod nvtx_utils.mod

vtaupsi_utils.mod.f90:$(SRCDIR)/vtaupsi_utils.mod.F90
vtaupsi_utils.mod.o:vtaupsi_utils.mod.f90 cppt.mod error_handling.mod \
                fftmain_utils.mod kinds.mod kpts.mod parac.mod \
                spin.mod system.mod tauf.mod tauofr_utils.mod \
                timer.mod zeroing_utils.mod

vtd2_utils.mod.f90:$(SRCDIR)/vtd2_utils.mod.F90
vtd2_utils.mod.o:vtd2_utils.mod.f90 dd_xc_utils.mod density_functionals_utils.mod \
                error_handling.mod func.mod kinds.mod linres.mod \
                lr_xcpot_utils.mod poin.mod spin.mod system.mod \
                timer.mod zeroing_utils.mod

wannier_center_utils.mod.f90:$(SRCDIR)/wannier_center_utils.mod.F90
wannier_center_utils.mod.o:wannier_center_utils.mod.f90 adat.mod \
                cdftmod.mod cnst.mod error_handling.mod fileopen_utils.mod \
                fileopenmod.mod hfxmod.mod ions.mod kinds.mod \
                meta_multiple_walkers_utils.mod metr.mod mm_dimmod.mod \
                mp_interface.mod mw.mod parac.mod pbc_utils.mod \
                pimd.mod response_pmod.mod ropt.mod store_types.mod \
                system.mod vdw_utils.mod vdwcmod.mod wann.mod

wannier_print_utils.mod.f90:$(SRCDIR)/wannier_print_utils.mod.F90
wannier_print_utils.mod.o:wannier_print_utils.mod.f90 cppt.mod \
                error_handling.mod fftmain_utils.mod kinds.mod \
                mp_interface.mod parac.mod readsr_utils.mod \
                system.mod wann.mod zeroing_utils.mod

wann.mod.f90:   $(SRCDIR)/wann.mod.F90
wann.mod.o:     wann.mod.f90 kinds.mod

wc_dos_utils.mod.f90:$(SRCDIR)/wc_dos_utils.mod.F90
wc_dos_utils.mod.o:wc_dos_utils.mod.f90 cnst.mod coor.mod error_handling.mod \
                fft_maxfft.mod fileopen_utils.mod fileopenmod.mod \
                forcep_utils.mod ions.mod kinds.mod meta_multiple_walkers_utils.mod \
                mp_interface.mod mw.mod ovlap_utils.mod parac.mod \
                prop.mod proylm_utils.mod spin.mod system.mod \
                timer.mod wann.mod zeroing_utils.mod fft_maxfft.mod

wfnio_utils.mod.f90:$(SRCDIR)/wfnio_utils.mod.F90
wfnio_utils.mod.o:wfnio_utils.mod.f90 bsym.mod error_handling.mod \
                geq0mod.mod gsortho_utils.mod io_utils.mod \
                kinds.mod kpclean_utils.mod kpts.mod mp_interface.mod \
                parac.mod prng_utils.mod randtowf_utils.mod \
                readmod.mod spin.mod system.mod timer.mod utils.mod \
                zeroing_utils.mod

wfn_print_utils.mod.f90:$(SRCDIR)/wfn_print_utils.mod.F90
wfn_print_utils.mod.o:wfn_print_utils.mod.f90 cppt.mod error_handling.mod \
                fileopen_utils.mod fileopenmod.mod g_loc.mod \
                kinds.mod mp_interface.mod parac.mod readsr_utils.mod \
                system.mod timer.mod zeroing_utils.mod kinds.mod \
                error_handling.mod timer.mod mp_interface.mod \
                system.mod parac.mod cppt.mod g_loc.mod fileopen_utils.mod \
                fileopenmod.mod readsr_utils.mod zeroing_utils.mod

wfopts_utils.mod.f90:$(SRCDIR)/wfopts_utils.mod.F90
wfopts_utils.mod.o:wfopts_utils.mod.f90 atwf.mod bswfo_utils.mod \
                cl_init_utils.mod clas.mod coor.mod cplngsmod.mod \
                ddip.mod do_perturbation_p_utils.mod elct.mod \
                error_handling.mod fnlalloc_utils.mod ions.mod \
                kinds.mod kpts.mod lr_in_utils.mod parac.mod \
                phfac_utils.mod prmem_utils.mod pslo.mod rwfopt_utils.mod \
                sfac.mod soft.mod spin.mod system.mod timer.mod \
                utils.mod vdwcmod.mod zeroing_utils.mod

wr30wfn_utils.mod.f90:$(SRCDIR)/wr30wfn_utils.mod.F90
wr30wfn_utils.mod.o:wr30wfn_utils.mod.f90 dotp_utils.mod error_handling.mod \
                kinds.mod summat_utils.mod timer.mod utils.mod \
                zeroing_utils.mod dotp_utils.mod kinds.mod \
                error_handling.mod timer.mod setbasis_utils.mod \
                mp_interface.mod system.mod parac.mod atwf.mod \
                ions.mod kpts.mod phfac_utils.mod ovlap_utils.mod \
                wr30wfn_utils.mod summat_utils.mod utils.mod \
                utils.mod zeroing_utils.mod dotp_utils.mod \
                kinds.mod error_handling.mod timer.mod mp_interface.mod \
                setbasis_utils.mod prng_utils.mod system.mod \
                parac.mod ions.mod atwf.mod kpts.mod spin.mod \
                gsortho_utils.mod phfac_utils.mod kinds.mod \
                error_handling.mod timer.mod mp_interface.mod \
                system.mod parac.mod kpts.mod io_utils.mod \
                zeroing_utils.mod cp_ieee_interface.mod kinds.mod \
                error_handling.mod timer.mod mp_interface.mod \
                system.mod parac.mod kpts.mod io_utils.mod \
                gsortho_utils.mod kinds.mod error_handling.mod \
                timer.mod mp_interface.mod system.mod parac.mod \
                kpts.mod bsym.mod utils.mod kinds.mod error_handling.mod \
                timer.mod mp_interface.mod prng_utils.mod system.mod \
                parac.mod kpts.mod spin.mod geq0mod.mod bsym.mod \
                io_utils.mod kpclean_utils.mod gsortho_utils.mod \
                randtowf_utils.mod zeroing_utils.mod utils.mod \
                kinds.mod error_handling.mod timer.mod kinds.mod \
                error_handling.mod timer.mod kinds.mod error_handling.mod \
                timer.mod mp_interface.mod system.mod parac.mod \
                geq0mod.mod rgs_utils.mod kinds.mod error_handling.mod \
                timer.mod mp_interface.mod system.mod parac.mod \
                part_1d.mod part_1d.mod io_utils.mod kinds.mod \
                error_handling.mod timer.mod mp_interface.mod \
                system.mod parac.mod part_1d.mod part_1d.mod \
                io_utils.mod zeroing_utils.mod

wrccfl_utils.mod.f90:$(SRCDIR)/wrccfl_utils.mod.F90
wrccfl_utils.mod.o:wrccfl_utils.mod.f90 bsym.mod fileopen_utils.mod \
                fileopenmod.mod kinds.mod parac.mod ropt.mod \
                system.mod

wrener_utils.mod.f90:$(SRCDIR)/wrener_utils.mod.F90
wrener_utils.mod.o:wrener_utils.mod.f90 adat.mod cdftmod.mod \
                clas.mod cnst.mod cnstpr_utils.mod conv.mod \
                coor.mod ener.mod fint.mod hubbardu.mod ions.mod \
                kinds.mod kpnt.mod kpts.mod machine.mod metr.mod \
                mm_dim_utils.mod mm_dimmod.mod mm_input.mod \
                parac.mod pslo.mod response_pmod.mod rmas.mod \
                ropt.mod spin.mod store_types.mod strs.mod \
                system.mod temps.mod timer.mod vdwcmod.mod \
                wrgeo_utils.mod

wrgeo_utils.mod.f90:$(SRCDIR)/wrgeo_utils.mod.F90
wrgeo_utils.mod.o:wrgeo_utils.mod.f90 adat.mod cnst.mod error_handling.mod \
                fileopen_utils.mod fileopenmod.mod ions.mod \
                kinds.mod mm_dim_utils.mod mm_dimmod.mod mm_input.mod \
                parac.mod ropt.mod store_types.mod system.mod \
                timer.mod

wrintf_utils.mod.f90:$(SRCDIR)/wrintf_utils.mod.F90
wrintf_utils.mod.o:wrintf_utils.mod.f90 error_handling.mod \
                fileopen_utils.mod fileopenmod.mod func.mod \
                ions.mod isos.mod kinds.mod kpnt.mod kpts.mod \
                parac.mod ragg.mod readsr_utils.mod store_types.mod \
                system.mod

write_pp_utils.mod.f90:$(SRCDIR)/write_pp_utils.mod.F90
write_pp_utils.mod.o:write_pp_utils.mod.f90 dpot.mod fileopen_utils.mod \
                fileopenmod.mod func.mod ions.mod kinds.mod \
                parac.mod sgpp.mod

wr_temps_utils.mod.f90:$(SRCDIR)/wr_temps_utils.mod.F90
wr_temps_utils.mod.o:wr_temps_utils.mod.f90 cnst.mod fileopen_utils.mod \
                fileopenmod.mod global_utils.mod kinds.mod \
                nose.mod parac.mod pimd.mod system.mod zeroing_utils.mod

wv30_utils.mod.f90:$(SRCDIR)/wv30_utils.mod.F90
wv30_utils.mod.o:wv30_utils.mod.f90 benc.mod cdft_utils.mod \
                cdftmod.mod cell.mod clas.mod cotr.mod elct.mod \
                ener.mod error_handling.mod fileopenmod.mod \
                filnmod.mod glemod.mod io_utils.mod ions.mod \
                kinds.mod kpnt.mod kpts.mod lscal.mod machine.mod \
                meta_multiple_walkers_utils.mod metr.mod mm_dimmod.mod \
                mm_extrap.mod mp_interface.mod mw.mod nose.mod \
                parac.mod pimd.mod poin.mod prng.mod readsr_utils.mod \
                rlbfgs_io.mod rw_linres_utils.mod shop_rest.mod \
                shop_rest_2.mod spin.mod store_types.mod string_utils.mod \
                symm.mod system.mod timer.mod wfnio_utils.mod \
                wrintf_utils.mod bicanonicalCpmd.mod kinds.mod \
                error_handling.mod timer.mod mp_interface.mod \
                system.mod parac.mod

xcener_utils.mod.f90:$(SRCDIR)/xcener_utils.mod.F90
xcener_utils.mod.o:xcener_utils.mod.f90 fft_maxfft.mod func.mod \
                functionals_utils.mod kinds.mod lsd_func_utils.mod \
                system.mod tbxc.mod timer.mod

xcstr_utils.mod.f90:$(SRCDIR)/xcstr_utils.mod.F90
xcstr_utils.mod.o:xcstr_utils.mod.f90 zeroing_utils.mod kinds.mod \
                error_handling.mod timer.mod system.mod parac.mod \
                cnst.mod cppt.mod strs.mod spin.mod geq0mod.mod \
                reshaper.mod fftnew_utils.mod fftmain_utils.mod \
                xcener_utils.mod zeroing_utils.mod

x_hjs.mod.f90:  $(SRCDIR)/x_hjs.mod.F90
x_hjs.mod.o:    x_hjs.mod.f90 error_handling.mod kinds.mod

xinr.mod.f90:   $(SRCDIR)/xinr.mod.F90
xinr.mod.o:     xinr.mod.f90 kinds.mod

ylmr2_utils.mod.f90:$(SRCDIR)/ylmr2_utils.mod.F90
ylmr2_utils.mod.o:ylmr2_utils.mod.f90 cnst.mod error_handling.mod \
                kinds.mod parac.mod

ylmr_utils.mod.f90:$(SRCDIR)/ylmr_utils.mod.F90
ylmr_utils.mod.o:ylmr_utils.mod.f90 bessm_utils.mod cnst.mod \
                error_handling.mod kinds.mod ylmr2_utils.mod

zdiis_utils.mod.f90:$(SRCDIR)/zdiis_utils.mod.F90
zdiis_utils.mod.o:zdiis_utils.mod.f90 kinds.mod linres.mod \
                mp_interface.mod odiis_utils.mod parac.mod \
                system.mod td_dav_utils.mod timer.mod zeroing_utils.mod

zeroing_utils.mod.f90:$(SRCDIR)/zeroing_utils.mod.F90
zeroing_utils.mod.o:zeroing_utils.mod.f90 azzero_utils.mod \
                kinds.mod

znum_mat_utils.mod.f90:$(SRCDIR)/znum_mat_utils.mod.F90
znum_mat_utils.mod.o:znum_mat_utils.mod.f90 cppt.mod error_handling.mod \
                geq0mod.mod kinds.mod metr.mod mp_interface.mod \
                parac.mod spin.mod system.mod zeroing_utils.mod

coordar.mod:    Gromos/coordar.mod.o
	@true
coordsz.mod:    Gromos/coordsz.mod.o
	@true
################################################################################
# Object dependencies: QM/MM modules
# 

Gromos/allocate_gromos.f:$(MODDIR)/Gromos/allocate_gromos.F
Gromos/allocate_gromos.o:Gromos/allocate_gromos.f kinds.mod \
                error_handling.mod prmem_utils.mod system.mod \
                parac.mod mm_parallel.mod coordsz.mod coordar.mod \
                $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h \
                $(MODDIR)/Gromos/forcesz.h $(MODDIR)/Gromos/nbpml.h \
                kinds.mod error_handling.mod coordsz.mod $(MODDIR)/Gromos/vector.h

Gromos/blockio.f:$(MODDIR)/Gromos/blockio.F
Gromos/blockio.o:Gromos/blockio.f $(MODDIR)/Gromos/dataid.h \
                kinds.mod error_handling.mod $(MODDIR)/Gromos/fileio.h \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                mp_interface.mod kinds.mod error_handling.mod \
                $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                mp_interface.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                mp_interface.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                mm_parallel.mod mp_interface.mod kinds.mod \
                error_handling.mod mm_parallel.mod mp_interface.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                mp_interface.mod kinds.mod error_handling.mod \
                mm_parallel.mod mp_interface.mod kinds.mod \
                error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                mp_interface.mod kinds.mod error_handling.mod \
                mm_parallel.mod mp_interface.mod coordsz.mod \
                kinds.mod error_handling.mod coordsz.mod kinds.mod \
                error_handling.mod kinds.mod error_handling.mod \
                mm_parallel.mod kinds.mod error_handling.mod \
                mm_parallel.mod kinds.mod error_handling.mod \
                mm_parallel.mod kinds.mod error_handling.mod \
                mm_parallel.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                mm_parallel.mod mp_interface.mod kinds.mod \
                error_handling.mod mm_parallel.mod kinds.mod \
                error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                mm_parallel.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                mm_parallel.mod mp_interface.mod kinds.mod \
                error_handling.mod kinds.mod error_handling.mod \
                $(MODDIR)/Gromos/disre.h $(MODDIR)/Gromos/formats.h \
                kinds.mod error_handling.mod mm_parallel.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                mm_parallel.mod mp_interface.mod kinds.mod \
                error_handling.mod kinds.mod error_handling.mod \
                mm_parallel.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                kinds.mod error_handling.mod coordsz.mod kinds.mod \
                error_handling.mod mm_parallel.mod coordsz.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                kinds.mod error_handling.mod mm_parallel.mod

Gromos/c06xxx.f:$(MODDIR)/Gromos/c06xxx.F
Gromos/c06xxx.o:Gromos/c06xxx.f kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod

Gromos/cobond.f:$(MODDIR)/Gromos/cobond.F
Gromos/cobond.o:Gromos/cobond.f kinds.mod error_handling.mod \
                mm_parallel.mod coordsz.mod $(MODDIR)/Gromos/box.h \
                $(MODDIR)/Gromos/pertsz.h $(MODDIR)/Gromos/cobond.h \
                kinds.mod error_handling.mod mm_parallel.mod \
                coordsz.mod kinds.mod error_handling.mod mm_parallel.mod \
                coordsz.mod kinds.mod error_handling.mod

Gromos/conat.f: $(MODDIR)/Gromos/conat.F
Gromos/conat.o: Gromos/conat.f kinds.mod error_handling.mod \
                coordsz.mod $(MODDIR)/Gromos/box.h

Gromos/coordar.mod.f90:$(MODDIR)/Gromos/coordar.mod.F90
Gromos/coordar.mod.o:Gromos/coordar.mod.f90 coordsz.mod

Gromos/coordsz.mod.f90:$(MODDIR)/Gromos/coordsz.mod.F90
Gromos/coordsz.mod.o:Gromos/coordsz.mod.f90

Gromos/fileio.f:$(MODDIR)/Gromos/fileio.F
Gromos/fileio.o:Gromos/fileio.f kinds.mod error_handling.mod \
                mm_parallel.mod mp_interface.mod $(MODDIR)/Gromos/units.h \
                kinds.mod error_handling.mod $(MODDIR)/Gromos/dataid.h \
                kinds.mod error_handling.mod mm_parallel.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                mm_parallel.mod error_handling.mod $(MODDIR)/Gromos/fileio.h \
                mm_parallel.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                mp_interface.mod kinds.mod error_handling.mod \
                mm_parallel.mod mp_interface.mod kinds.mod \
                error_handling.mod mm_parallel.mod kinds.mod \
                error_handling.mod mm_parallel.mod mp_interface.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                mm_parallel.mod kinds.mod error_handling.mod \
                mm_parallel.mod kinds.mod error_handling.mod \
                mm_parallel.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod parac.mod mp_interface.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod

Gromos/gauss.f: $(MODDIR)/Gromos/gauss.F
Gromos/gauss.o: Gromos/gauss.f kinds.mod error_handling.mod

Gromos/latsum.f:$(MODDIR)/Gromos/latsum.F
Gromos/latsum.o:Gromos/latsum.f kinds.mod error_handling.mod \
                mm_parallel.mod kinds.mod error_handling.mod \
                mm_parallel.mod kinds.mod error_handling.mod \
                mm_parallel.mod kinds.mod error_handling.mod \
                mm_parallel.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod

Gromos/mdutils.f:$(MODDIR)/Gromos/mdutils.F
Gromos/mdutils.o:Gromos/mdutils.f kinds.mod error_handling.mod \
                mm_parallel.mod parac.mod mp_interface.mod \
                coordsz.mod coordar.mod $(MODDIR)/Gromos/toposz.h \
                $(MODDIR)/Gromos/units.h $(MODDIR)/Gromos/topoar.h \
                $(MODDIR)/Gromos/md.h $(MODDIR)/Gromos/box.h \
                kinds.mod error_handling.mod coordsz.mod $(MODDIR)/Gromos/cobond.h \
                kinds.mod error_handling.mod coordsz.mod coordar.mod \
                kinds.mod error_handling.mod $(MODDIR)/Gromos/pertsz.h \
                $(MODDIR)/Gromos/pertar.h kinds.mod error_handling.mod \
                coordsz.mod $(MODDIR)/Gromos/forcesz.h kinds.mod \
                error_handling.mod kinds.mod error_handling.mod \
                coordsz.mod coordar.mod $(MODDIR)/Gromos/forcear.h \
                kinds.mod error_handling.mod coordsz.mod coordar.mod \
                kinds.mod error_handling.mod coordar.mod kinds.mod \
                error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod

Gromos/mm_concg.f:$(MODDIR)/Gromos/mm_concg.F
Gromos/mm_concg.o:Gromos/mm_concg.f kinds.mod error_handling.mod \
                coordsz.mod $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h \
                $(MODDIR)/Gromos/box.h kinds.mod error_handling.mod \
                coordsz.mod

Gromos/mm_force.f:$(MODDIR)/Gromos/mm_force.F
Gromos/mm_force.o:Gromos/mm_force.f kinds.mod error_handling.mod \
                ropt.mod mm_parallel.mod mp_interface.mod azzero_utils.mod \
                coordar.mod $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h \
                $(MODDIR)/Gromos/md.h $(MODDIR)/Gromos/box.h \
                $(MODDIR)/Gromos/formats.h $(MODDIR)/Gromos/forcesz.h \
                $(MODDIR)/Gromos/forcear.h $(MODDIR)/Gromos/runmd.h \
                $(MODDIR)/Gromos/pertsz.h $(MODDIR)/Gromos/pertar.h \
                $(MODDIR)/Gromos/mm_save_config.h

Gromos/mm_for.f:$(MODDIR)/Gromos/mm_for.F
Gromos/mm_for.o:Gromos/mm_for.f kinds.mod error_handling.mod \
                mm_input.mod mm_parallel.mod coordsz.mod $(MODDIR)/Gromos/md.h \
                $(MODDIR)/Gromos/box.h $(MODDIR)/Gromos/toposz.h \
                $(MODDIR)/Gromos/topoar.h $(MODDIR)/Gromos/forcesz.h \
                $(MODDIR)/Gromos/forcear.h $(MODDIR)/Gromos/cobond.h \
                $(MODDIR)/Gromos/restx.h $(MODDIR)/Gromos/pertsz.h \
                $(MODDIR)/Gromos/pertar.h $(MODDIR)/Gromos/disre.h \
                kinds.mod error_handling.mod coordsz.mod kinds.mod \
                error_handling.mod coordsz.mod kinds.mod error_handling.mod \
                coordsz.mod mm_input.mod kinds.mod error_handling.mod \
                mm_input.mod coordsz.mod kinds.mod error_handling.mod \
                coordsz.mod

Gromos/mm_restx.f:$(MODDIR)/Gromos/mm_restx.F
Gromos/mm_restx.o:Gromos/mm_restx.f kinds.mod error_handling.mod \
                mm_parallel.mod mp_interface.mod coordsz.mod \
                $(MODDIR)/Gromos/box.h $(MODDIR)/Gromos/restx.h

Gromos/mm_setup_dr.f:$(MODDIR)/Gromos/mm_setup_dr.F
Gromos/mm_setup_dr.o:Gromos/mm_setup_dr.f kinds.mod error_handling.mod \
                mm_parallel.mod coordar.mod $(MODDIR)/Gromos/toposz.h \
                $(MODDIR)/Gromos/topoar.h $(MODDIR)/Gromos/md.h \
                $(MODDIR)/Gromos/box.h $(MODDIR)/Gromos/formats.h \
                $(MODDIR)/Gromos/forcesz.h $(MODDIR)/Gromos/forcear.h \
                $(MODDIR)/Gromos/runmd.h $(MODDIR)/Gromos/pertsz.h \
                $(MODDIR)/Gromos/pertar.h $(MODDIR)/Gromos/mm_save_config.h \
                kinds.mod error_handling.mod coordsz.mod kinds.mod \
                error_handling.mod coordsz.mod coordsz.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                coordsz.mod kinds.mod error_handling.mod coordsz.mod \
                kinds.mod error_handling.mod coordsz.mod kinds.mod \
                error_handling.mod coordsz.mod kinds.mod error_handling.mod \
                coordsz.mod kinds.mod error_handling.mod kinds.mod \
                error_handling.mod kinds.mod error_handling.mod \
                coordsz.mod kinds.mod error_handling.mod coordsz.mod \
                kinds.mod error_handling.mod coordar.mod kinds.mod \
                error_handling.mod mm_parallel.mod coordsz.mod

Gromos/mm_setup.f:$(MODDIR)/Gromos/mm_setup.F
Gromos/mm_setup.o:Gromos/mm_setup.f kinds.mod error_handling.mod \
                system.mod parac.mod mm_input.mod fileopenmod.mod \
                fileopen_utils.mod mm_parallel.mod mp_interface.mod \
                coordar.mod $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h \
                $(MODDIR)/Gromos/forcesz.h $(MODDIR)/Gromos/box.h \
                $(MODDIR)/Gromos/md.h $(MODDIR)/Gromos/units.h \
                kinds.mod error_handling.mod mm_parallel.mod \
                coordsz.mod

Gromos/mm_shake.f:$(MODDIR)/Gromos/mm_shake.F
Gromos/mm_shake.o:Gromos/mm_shake.f kinds.mod error_handling.mod \
                coordsz.mod $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/box.h \
                kinds.mod error_handling.mod coordsz.mod

Gromos/nbpml.f: $(MODDIR)/Gromos/nbpml.F
Gromos/nbpml.o: Gromos/nbpml.f kinds.mod error_handling.mod \
                mm_parallel.mod coordsz.mod $(MODDIR)/Gromos/toposz.h \
                $(MODDIR)/Gromos/topoar.h $(MODDIR)/Gromos/box.h \
                $(MODDIR)/Gromos/forcesz.h $(MODDIR)/Gromos/pertsz.h \
                $(MODDIR)/Gromos/pertar.h $(MODDIR)/Gromos/nbpml.h \
                kinds.mod error_handling.mod coordsz.mod kinds.mod \
                error_handling.mod mm_parallel.mod mm_input.mod \
                mp_interface.mod coordsz.mod kinds.mod error_handling.mod \
                mm_parallel.mod coordsz.mod kinds.mod error_handling.mod \
                mm_parallel.mod coordsz.mod $(MODDIR)/Gromos/vector.h \
                kinds.mod error_handling.mod mm_parallel.mod \
                coordsz.mod kinds.mod error_handling.mod mm_parallel.mod \
                coordsz.mod kinds.mod error_handling.mod coordsz.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod

Gromos/nonbml.f:$(MODDIR)/Gromos/nonbml.F
Gromos/nonbml.o:Gromos/nonbml.f kinds.mod error_handling.mod \
                mm_input.mod mm_parallel.mod coordsz.mod $(MODDIR)/Gromos/box.h \
                $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h \
                $(MODDIR)/Gromos/pertsz.h $(MODDIR)/Gromos/pertar.h \
                $(MODDIR)/Gromos/forcesz.h kinds.mod error_handling.mod \
                mm_input.mod mm_parallel.mod coordsz.mod $(MODDIR)/Gromos/vector.h

Gromos/posio.f: $(MODDIR)/Gromos/posio.F
Gromos/posio.o: Gromos/posio.f kinds.mod error_handling.mod \
                coordsz.mod kinds.mod error_handling.mod mm_parallel.mod \
                coordsz.mod $(MODDIR)/Gromos/dataid.h kinds.mod \
                error_handling.mod mm_parallel.mod mp_interface.mod \
                coordar.mod $(MODDIR)/Gromos/fileio.h kinds.mod \
                error_handling.mod mm_parallel.mod kinds.mod \
                error_handling.mod coordar.mod kinds.mod error_handling.mod \
                coordar.mod

Gromos/random.f:$(MODDIR)/Gromos/random.F
Gromos/random.o:Gromos/random.f kinds.mod error_handling.mod

Gromos/rdmd.f:  $(MODDIR)/Gromos/rdmd.F
Gromos/rdmd.o:  Gromos/rdmd.f kinds.mod error_handling.mod \
                mm_parallel.mod parac.mod mp_interface.mod \
                $(MODDIR)/Gromos/md.h $(MODDIR)/Gromos/mdblock.h \
                $(MODDIR)/Gromos/fileio.h kinds.mod error_handling.mod \
                coordar.mod $(MODDIR)/Gromos/box.h $(MODDIR)/Gromos/formats.h \
                $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h \
                $(MODDIR)/Gromos/forcesz.h kinds.mod error_handling.mod \
                coordsz.mod kinds.mod error_handling.mod kinds.mod \
                error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                mm_parallel.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                coordar.mod kinds.mod error_handling.mod coordar.mod \
                kinds.mod error_handling.mod parac.mod mp_interface.mod \
                coordar.mod kinds.mod error_handling.mod coordar.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod coordsz.mod kinds.mod \
                error_handling.mod coordsz.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod coordsz.mod kinds.mod \
                error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                mm_parallel.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod coordsz.mod kinds.mod \
                error_handling.mod coordsz.mod kinds.mod error_handling.mod \
                coordsz.mod kinds.mod error_handling.mod kinds.mod \
                error_handling.mod $(MODDIR)/Gromos/inc_rdmd5.h \
                kinds.mod error_handling.mod $(MODDIR)/Gromos/inc_rdmd4.h

Gromos/rdtopo.f:$(MODDIR)/Gromos/rdtopo.F
Gromos/rdtopo.o:Gromos/rdtopo.f kinds.mod error_handling.mod \
                $(MODDIR)/Gromos/units.h kinds.mod error_handling.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                mp_interface.mod $(MODDIR)/Gromos/toposz.h \
                $(MODDIR)/Gromos/topoar.h $(MODDIR)/Gromos/fileio.h \
                $(MODDIR)/Gromos/topblock.h kinds.mod error_handling.mod \
                mm_parallel.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                mp_interface.mod $(MODDIR)/Gromos/formats.h \
                kinds.mod error_handling.mod mm_parallel.mod \
                mp_interface.mod kinds.mod error_handling.mod \
                mm_parallel.mod mp_interface.mod kinds.mod \
                error_handling.mod mm_parallel.mod mp_interface.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                mp_interface.mod kinds.mod error_handling.mod \
                mm_parallel.mod mp_interface.mod kinds.mod \
                error_handling.mod mm_parallel.mod mp_interface.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                mp_interface.mod kinds.mod error_handling.mod \
                mm_parallel.mod mp_interface.mod kinds.mod \
                error_handling.mod mm_parallel.mod mp_interface.mod \
                kinds.mod error_handling.mod mm_parallel.mod \
                mp_interface.mod kinds.mod error_handling.mod \
                mm_parallel.mod mp_interface.mod kinds.mod \
                error_handling.mod mm_parallel.mod mp_interface.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                mm_parallel.mod

Gromos/schifezza.f:$(MODDIR)/Gromos/schifezza.F
Gromos/schifezza.o:Gromos/schifezza.f kinds.mod error_handling.mod \
                coordsz.mod coordar.mod $(MODDIR)/Gromos/box.h \
                $(MODDIR)/Gromos/dataid.h $(MODDIR)/Gromos/fileio.h \
                $(MODDIR)/Gromos/forcesz.h $(MODDIR)/Gromos/formats.h \
                $(MODDIR)/Gromos/md.h $(MODDIR)/Gromos/mdblock.h \
                $(MODDIR)/Gromos/pertsz.h $(MODDIR)/Gromos/pertar.h \
                $(MODDIR)/Gromos/ptblock.h $(MODDIR)/Gromos/toposz.h \
                $(MODDIR)/Gromos/topoar.h $(MODDIR)/Gromos/units.h \
                kinds.mod error_handling.mod coordsz.mod coordar.mod \
                $(MODDIR)/Gromos/topblock.h

Gromos/string.f:$(MODDIR)/Gromos/string.F
Gromos/string.o:Gromos/string.f kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod

Gromos/traco.f: $(MODDIR)/Gromos/traco.F
Gromos/traco.o: Gromos/traco.f kinds.mod error_handling.mod \
                coordsz.mod $(MODDIR)/Gromos/box.h

Gromos/wrtopo.f:$(MODDIR)/Gromos/wrtopo.F
Gromos/wrtopo.o:Gromos/wrtopo.f kinds.mod error_handling.mod \
                $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h \
                $(MODDIR)/Gromos/topblock.h kinds.mod error_handling.mod \
                kinds.mod error_handling.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod

MM_Interface/mm_add_dummy.f:$(MODDIR)/MM_Interface/mm_add_dummy.F
MM_Interface/mm_add_dummy.o:MM_Interface/mm_add_dummy.f kinds.mod \
                error_handling.mod mm_input.mod coordar.mod \
                $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h \
                $(MODDIR)/Gromos/md.h $(MODDIR)/Gromos/forcesz.h \
                $(MODDIR)/Gromos/forcear.h

MM_Interface/mm_add_hydrogen.f:$(MODDIR)/MM_Interface/mm_add_hydrogen.F
MM_Interface/mm_add_hydrogen.o:MM_Interface/mm_add_hydrogen.f \
                kinds.mod error_handling.mod mm_input.mod coordar.mod \
                $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h \
                kinds.mod error_handling.mod constr_utils.mod

MM_Interface/mm_cap_H.f:$(MODDIR)/MM_Interface/mm_cap_H.F
MM_Interface/mm_cap_H.o:MM_Interface/mm_cap_H.f kinds.mod error_handling.mod \
                mm_input.mod system.mod mm_parallel.mod coordar.mod \
                $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h

MM_Interface/mm_center.f:$(MODDIR)/MM_Interface/mm_center.F
MM_Interface/mm_center.o:MM_Interface/mm_center.f kinds.mod \
                error_handling.mod system.mod parac.mod cell.mod \
                rmas.mod mm_dimmod.mod isos.mod mm_input.mod \
                bsym.mod

MM_Interface/mm_charge_rest.f:$(MODDIR)/MM_Interface/mm_charge_rest.F
MM_Interface/mm_charge_rest.o:MM_Interface/mm_charge_rest.f \
                kinds.mod error_handling.mod system.mod parac.mod \
                mm_dimmod.mod mm_input.mod

MM_Interface/mm_charges.f:$(MODDIR)/MM_Interface/mm_charges.F
MM_Interface/mm_charges.o:MM_Interface/mm_charges.f kinds.mod \
                error_handling.mod mp_interface.mod prmem_utils.mod \
                system.mod parac.mod cell.mod mm_input.mod \
                atwf.mod mm_dimmod.mod ions.mod azzero_utils.mod

MM_Interface/mm_detit.f:$(MODDIR)/MM_Interface/mm_detit.F
MM_Interface/mm_detit.o:MM_Interface/mm_detit.f kinds.mod error_handling.mod \
                system.mod parac.mod ions.mod mm_input.mod \
                fileopenmod.mod fileopen_utils.mod cotr.mod \
                coordsz.mod $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h

MM_Interface/mm_elstat_atdens.f:$(MODDIR)/MM_Interface/mm_elstat_atdens.F
MM_Interface/mm_elstat_atdens.o:MM_Interface/mm_elstat_atdens.f \
                kinds.mod error_handling.mod mp_interface.mod \
                system.mod parac.mod mm_dimmod.mod ions.mod \
                azzero_utils.mod kinds.mod error_handling.mod \
                system.mod atwf.mod mm_dimmod.mod kinds.mod \
                error_handling.mod kinds.mod error_handling.mod

MM_Interface/mm_elstat.f:$(MODDIR)/MM_Interface/mm_elstat.F
MM_Interface/mm_elstat.o:MM_Interface/mm_elstat.f kinds.mod \
                error_handling.mod system.mod parac.mod machine.mod \
                mp_interface.mod prmem_utils.mod epot_types.mod \
                cnst.mod cell.mod adat.mod ions.mod ropt.mod \
                store_types.mod azzero_utils.mod mm_dimmod.mod \
                mm_input.mod mm_parallel.mod fileopenmod.mod \
                fileopen_utils.mod forcematch.mod efld.mod

MM_Interface/mm_elstat_sr.f:$(MODDIR)/MM_Interface/mm_elstat_sr.F
MM_Interface/mm_elstat_sr.o:MM_Interface/mm_elstat_sr.f kinds.mod \
                error_handling.mod timer.mod mp_interface.mod \
                system.mod parac.mod epot_types.mod cell.mod \
                mm_input.mod azzero_utils.mod

MM_Interface/mm_excl.f:$(MODDIR)/MM_Interface/mm_excl.F
MM_Interface/mm_excl.o:MM_Interface/mm_excl.f kinds.mod error_handling.mod \
                timer.mod system.mod parac.mod mm_dimmod.mod \
                mm_input.mod cell.mod azzero_utils.mod

MM_Interface/mm_excl_mech.f:$(MODDIR)/MM_Interface/mm_excl_mech.F
MM_Interface/mm_excl_mech.o:MM_Interface/mm_excl_mech.f kinds.mod \
                error_handling.mod timer.mod system.mod parac.mod \
                mm_dimmod.mod mm_input.mod azzero_utils.mod

MM_Interface/mm_flex_solv.f:$(MODDIR)/MM_Interface/mm_flex_solv.F
MM_Interface/mm_flex_solv.o:MM_Interface/mm_flex_solv.f kinds.mod \
                error_handling.mod system.mod mm_input.mod \
                store_types.mod fileopenmod.mod fileopen_utils.mod \
                coordar.mod $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h \
                $(MODDIR)/Gromos/md.h

MM_Interface/mm_gcp.f:$(MODDIR)/MM_Interface/mm_gcp.F
MM_Interface/mm_gcp.o:MM_Interface/mm_gcp.f kinds.mod error_handling.mod \
                mm_dimmod.mod cnst.mod coordsz.mod kinds.mod \
                error_handling.mod mm_dimmod.mod cnst.mod coordsz.mod \
                kinds.mod error_handling.mod cnst.mod mm_dimmod.mod \
                coordsz.mod

MM_Interface/mm_get_NSX.f:$(MODDIR)/MM_Interface/mm_get_NSX.F
MM_Interface/mm_get_NSX.o:MM_Interface/mm_get_NSX.f kinds.mod \
                error_handling.mod system.mod ions.mod mm_input.mod \
                mm_parallel.mod coordsz.mod $(MODDIR)/Gromos/toposz.h \
                $(MODDIR)/Gromos/topoar.h

MM_Interface/mm_invert.f:$(MODDIR)/MM_Interface/mm_invert.F
MM_Interface/mm_invert.o:MM_Interface/mm_invert.f kinds.mod \
                error_handling.mod timer.mod

MM_Interface/mm_long_range_classic.f:$(MODDIR)/MM_Interface/mm_long_range_classic.F
MM_Interface/mm_long_range_classic.o:MM_Interface/mm_long_range_classic.f \
                kinds.mod error_handling.mod system.mod parac.mod \
                mm_dimmod.mod ions.mod azzero_utils.mod

MM_Interface/mm_long_range_potential.f:$(MODDIR)/MM_Interface/mm_long_range_potential.F
MM_Interface/mm_long_range_potential.o:MM_Interface/mm_long_range_potential.f \
                kinds.mod error_handling.mod system.mod parac.mod \
                cell.mod

MM_Interface/mm_min_im.f:$(MODDIR)/MM_Interface/mm_min_im.F
MM_Interface/mm_min_im.o:MM_Interface/mm_min_im.f kinds.mod \
                error_handling.mod mm_dimmod.mod cell.mod

MM_Interface/mm_multipole.f:$(MODDIR)/MM_Interface/mm_multipole.F
MM_Interface/mm_multipole.o:MM_Interface/mm_multipole.f kinds.mod \
                error_handling.mod mp_interface.mod system.mod \
                parac.mod cell.mod mm_dimmod.mod ions.mod

MM_Interface/mm_nlist.f:$(MODDIR)/MM_Interface/mm_nlist.F
MM_Interface/mm_nlist.o:MM_Interface/mm_nlist.f kinds.mod error_handling.mod \
                timer.mod system.mod parac.mod mm_dimmod.mod \
                ions.mod cell.mod adat.mod mm_input.mod fileopenmod.mod \
                fileopen_utils.mod forcematch.mod

MM_Interface/mm_pardef.f:$(MODDIR)/MM_Interface/mm_pardef.F
MM_Interface/mm_pardef.o:MM_Interface/mm_pardef.f kinds.mod \
                error_handling.mod adat.mod cnst.mod coordsz.mod \
                $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h \
                $(MODDIR)/Gromos/box.h

MM_Interface/mm_quantum_topo.f:$(MODDIR)/MM_Interface/mm_quantum_topo.F
MM_Interface/mm_quantum_topo.o:MM_Interface/mm_quantum_topo.f \
                kinds.mod error_handling.mod prmem_utils.mod \
                system.mod parac.mod mm_input.mod forcematch.mod \
                mm_parallel.mod forcematch_utils.mod forcematch_kfit_utils.mod \
                $(MODDIR)/Gromos/toposz.h $(MODDIR)/Gromos/topoar.h

MM_Interface/mm_readgromos.f:$(MODDIR)/MM_Interface/mm_readgromos.F
MM_Interface/mm_readgromos.o:MM_Interface/mm_readgromos.f kinds.mod \
                error_handling.mod coordar.mod kinds.mod error_handling.mod \
                coordar.mod

MM_Interface/mm_read_qmmm_input.f:$(MODDIR)/MM_Interface/mm_read_qmmm_input.F
MM_Interface/mm_read_qmmm_input.o:MM_Interface/mm_read_qmmm_input.f \
                kinds.mod error_handling.mod mp_interface.mod \
                system.mod parac.mod mm_input.mod mm_dimmod.mod \
                efld.mod forcematch.mod readsr_utils.mod inscan_utils.mod \
                coordsz.mod $(MODDIR)/Gromos/vector.h $(MODDIR)/Gromos/forcesz.h \
                $(MODDIR)/Gromos/toposz.h kinds.mod error_handling.mod \
                system.mod readsr_utils.mod coordsz.mod $(MODDIR)/Gromos/box.h \
                kinds.mod error_handling.mod mm_input.mod forcematch.mod

MM_Interface/mm_reflex.f:$(MODDIR)/MM_Interface/mm_reflex.F
MM_Interface/mm_reflex.o:MM_Interface/mm_reflex.f kinds.mod \
                error_handling.mod system.mod ions.mod mm_dimmod.mod \
                isos.mod cell.mod mm_input.mod comvel_utils.mod

MM_Interface/mm_restart_traj.f:$(MODDIR)/MM_Interface/mm_restart_traj.F
MM_Interface/mm_restart_traj.o:MM_Interface/mm_restart_traj.f \
                kinds.mod error_handling.mod system.mod mp_interface.mod \
                parac.mod ions.mod store_types.mod mm_input.mod \
                mm_dim_utils.mod mm_dimmod.mod fileopenmod.mod \
                fileopen_utils.mod

MM_Interface/mm_short_range_classic.f:$(MODDIR)/MM_Interface/mm_short_range_classic.F
MM_Interface/mm_short_range_classic.o:MM_Interface/mm_short_range_classic.f \
                kinds.mod error_handling.mod system.mod mm_dimmod.mod \
                ions.mod mm_ion_dens.mod kinds.mod error_handling.mod

MM_Interface/mm_solv_constr.f:$(MODDIR)/MM_Interface/mm_solv_constr.F
MM_Interface/mm_solv_constr.o:MM_Interface/mm_solv_constr.f \
                kinds.mod error_handling.mod prmem_utils.mod \
                system.mod rmas.mod mm_dimmod.mod cnst.mod \
                azzero_utils.mod kinds.mod error_handling.mod \
                kinds.mod error_handling.mod

MM_Interface/mm_transl_c0.f:$(MODDIR)/MM_Interface/mm_transl_c0.F
MM_Interface/mm_transl_c0.o:MM_Interface/mm_transl_c0.f kinds.mod \
                error_handling.mod system.mod parac.mod cppt.mod \
                gvec.mod isos.mod mm_dimmod.mod geq0mod.mod \
                elct.mod mm_extrap.mod

MM_Interface/mm_write_gromos_coord.f:$(MODDIR)/MM_Interface/mm_write_gromos_coord.F
MM_Interface/mm_write_gromos_coord.o:MM_Interface/mm_write_gromos_coord.f \
                kinds.mod error_handling.mod fileopenmod.mod \
                fileopen_utils.mod machine.mod coordar.mod \
                $(MODDIR)/Gromos/md.h $(MODDIR)/Gromos/toposz.h \
                $(MODDIR)/Gromos/topoar.h $(MODDIR)/Gromos/box.h

MM_Interface/mm_write_potential.f:$(MODDIR)/MM_Interface/mm_write_potential.F
MM_Interface/mm_write_potential.o:MM_Interface/mm_write_potential.f \
                kinds.mod error_handling.mod mp_interface.mod \
                system.mod parac.mod cell.mod fileopenmod.mod \
                fileopen_utils.mod


# create a regular timestamped cpmd package from the cvs tree
gtar-cpmd  :
	@( d=`date +%Y%m%d` ; b=`basename $$PWD` ; cd ${CPMDROOT}/.. ;	\
	tar -c --exclude \*,v --exclude \*.bak --exclude \*~ --exclude .svn	\
	--exclude \*.o --exclude \*.a --exclude \*.x --exclude \*.log 		\
	--exclude \*.out --exclude \*.prj --exclude \*.chk			\
	--exclude modules/MM_Interface --exclude modules/Gromos 		\
        --exclude modules/IPhigenie_Interface					\
	--exclude \*.orig --exclude \*.rej -zvvf $$b-$$d.tar.gz $$b &&	\
	echo successfully created $$b-$$d-tar.gz ; cd $$b )

# create a QM/MM timestamped cpmd package from the cvs tree
gtar-qmmm :
	@( d=`date +%Y%m%d` ; b=`basename $$PWD` ; cd ${CPMDROOT}/.. ;	\
	tar -c --exclude \*,v --exclude \*.bak --exclude \*~ --exclude .svn	\
	--exclude \*.o --exclude \*.a --exclude \*.x --exclude \*.log 		\
	--exclude \*.out --exclude \*.prj --exclude \*.chk			\
	--exclude \*.orig --exclude \*.rej -zvvf $$b-$$d.tar.gz $$b &&	\
	echo successfully created $$b-$$d-tar.gz ; cd $$b )

