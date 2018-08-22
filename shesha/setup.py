from setuptools import setup

setup(
        name='shesha',
        version='3.0.0',
        # author=['Arnaud Sevin'],
        # author_email=['arnaud.sevin@obspm.fr'],
        # description='',
        # long_description='',
        packages=[
                'data', 'data.par', 'data.par.par4bench', 'shesha', 'shesha.ao',
                'shesha.config', 'shesha.init', 'shesha.scripts', 'shesha.sim',
                'shesha.supervisor', 'shesha.util', 'shesha.widgets'
        ],
        # packages=find_packages("shesha"),
        package_dir={
                'data': 'shesha/data',
                # 'data.par': 'shesha/data/par',
                # 'data.par.par4bench': 'shesha/data/par/par4bench',
                'shesha': 'shesha/shesha',
                # 'shesha.ao': 'shesha/shesha/ao',
                # 'shesha.config': 'shesha/shesha/config',
                # 'shesha.init': 'shesha/shesha/init',
                # 'shesha.scripts': 'shesha/shesha/scripts',
                # 'shesha.sim': 'shesha/shesha/sim',
                # 'shesha.supervisor': 'shesha/shesha/supervisor',
                # 'shesha.sutra_pybind': 'shesha/shesha/sutra_pybind',
                # 'shesha.util': 'shesha/shesha/util',
                # 'shesha.widgets': 'shesha/shesha/widgets',
        },
        package_data={'data': ['layouts/SCAO_PYR.area', 'layouts/SCAO_SH.area']},
        include_package_data=True,
        # install_requires=['compass_sim'],
        zip_safe=False, )
