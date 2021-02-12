
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

data_files = ['README.md', 'LICENSE', '*.ipynb']

setuptools.setup(
    include_package_data=True,
    name="wpspec",
    #packages=['wpspec'],
    packages=setuptools.find_packages(),
    package_data={'wpspec':data_files},
    package_dir={'wpspec': 'wpspec'},
    #packages=setuptools.find_packages(),
    #package_dir={'wpspec': 'wpspec'},
    #package_data={'wpspec': 'wpspec'},
    version="0.01.b",
    author="Foley Lab",
    author_email="foleyj10@wpunj.edu",
    description="A pedagogical package for simulating light and matter.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://foleylab.github.io/",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: MacOS :: MacOS X",
	"Operating System :: Microsoft :: Windows :: Windows 10",
	"Operating System :: POSIX :: Linux",
    ],
)
