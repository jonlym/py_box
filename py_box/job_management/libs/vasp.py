import os


def setup_vasp_parser(parser):
    parser.add_argument('-g', '--gamma', dest='gamma', default=False,
                        action='store_true', 
                        help='Turns on gamma point version of VASP')
    parser.add_argument('-l', '--noncollinear', dest='collinear', default=False,
                        action='store_true', help='Uses the non-collinear version of VASP')

    help = 'Which version of the VASP PAW\'s to use:  '
    help += 'old, v52, v54 (default)'
    parser.add_argument('-p', '--paw', dest='paw', default='v54',
                        type=str, help=help)

    help = 'Which version of VASP to use:  '
    help += '5.3, 5.4, 5.4_sol, 5.4_beef, 5.4_nbo'
    parser.add_argument('-v', '--version', dest='version', default='5.4',
                        type=str, help=help)

    help = 'Controls whether or not a node is exclusive to you.\n'
    parser.add_argument('-x', '--exclusive', dest='exclusive', default=True,
                        action='store_true', help=help)

def setup_vasp_environment(args, environment):
#    print 'Setting up VASP environment'
    project_home = os.environ['WORKDIR']
    shell = os.environ['SHELL'].split('/')[-1].lower()

    if args.version.lower() == '5.3':
        environment += 'vpkg_require vasp/5.3.2:d3,gamma_only,intel,mpi,vtst python-numpy python-scipy,vpkg_require python-pandas/python2.7.8\n\n'

        if args.gamma:
            print 'Using gamma-point version of VASP'
            exe = 'mpiexec -n %i %s/sw/vasp/5.3.2-intel64-openmpi+VTST+D3+GAMMA/vasp' % (args.nproc, project_home)
        elif args.collinear:
            print 'Using non-collinear version of VASP'
            exe = 'mpiexec -n %i %s/programs/vasp/vasp.5.2-no-collinear/vasp' % (args.nproc, project_home)
        elif not args.gamma and not args.collinear:
            print 'Using standard version of VASP'
            exe = 'mpiexec -n %i %s/sw/vasp/5.3.2-intel64-openmpi+VTST+D3/vasp' % (args.nproc, project_home)
    elif args.version.lower() in ['5.4', '5.4_sol', '5.4_beef']:
        if '_' in args.version:
            version = args.version.split('_')[1]
            vasp_base = '%s/programs/vasp5.4/vasp.5.4.1_%s' % (project_home, version)
        else:
            version = 'standard'
            vasp_base = '%s/programs/vasp5.4/vasp.5.4.1' % (project_home)

        environment += 'vpkg_require openmpi/1.10.2-intel64-2016 python-numpy python-scipy\n\n'
        vasp_exe_base = 'mpiexec -n %i %s' % (args.nproc, vasp_base)

        if args.gamma:
            print 'Using gamma-point version of VASP5.4 (%s)' % version
            exe = '%s/build/gam/vasp' % vasp_exe_base 
        elif args.collinear:
            print 'Using non-collinear version of VASP5.4 (%s)' % version
            exe = '%s/build/ncl/vasp' % vasp_exe_base
        elif not args.gamma and not args.collinear:
            print 'Using standard version of VASP5.4 (%s)' % version
            exe = '%s/build/std/vasp' % vasp_exe_base

    if args.paw == 'old':
        vasp_pp_path = '%s/programs/vasp_psp/old/\n' % project_home
    else:
        vasp_pp_path = '%s/programs/vasp_psp/%s/\n' % (project_home, args.paw)

    if 'csh' in shell:
        environment += 'setenv VASP_COMMAND \"%s\"\n' % exe
        environment += 'setenv VASP_PP_PATH %s\n' % vasp_pp_path
    elif 'bash' in shell:
        environment += 'export VASP_COMMAND=\"%s\"\n' % exe
        environment += 'export VASP_PP_PATH=%s\n' % vasp_pp_path

    if args.nproc % 20 > 0:
        msg = 'VASP requires all the cores on a node to be used.\nPlease set the number of cores '
        msg += 'to be used to a multiple of 20.'
        raise RuntimeError(msg) 

    return environment
