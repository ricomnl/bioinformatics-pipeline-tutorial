"""setup.py"""
from setuptools import setup

setup(
    name="bioinformatics_pipeline_tutorial",
    version="0.0.1",
    author="Rico Meinl",
    author_email="dev@rmeinl.com",
    description="This is the accompanying GitHub repository for this blog post: https://ricomnl.com/blog/bottom-up-bioinformatics-pipeline.",
    license="MIT",
    keywords="bioinformatics tutorial",
    url="https://github.com/ricomnl/bioinformatics-pipeline-tutorial",
    packages=["bioinformatics_pipeline_tutorial"],
    install_requires=["redun", "plotly", "kaleido"],
)
