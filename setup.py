from setuptools import setup
import glob

setup(name='pipe',
      version='0.0.1',
      packages=['pipe','pipe.slurmpy', 'pipe.utils'],
      description='NGS Pipeline',
      url='https://github.com/AndersenLab/vcf-toolbox',
      author='Daniel Cook',
      author_email='danielecook@gmail.com',
      license='MIT',
      entry_points="""
      [console_scripts]
      pipe= pipe.pipe:main
      """,
      install_requires=["cython","docopt", "cyvcf2", "biopython", "clint", "requests"],
      data_files=[glob.glob("scripts/*")],
      zip_safe=False)