from setuptools import setup


setup(
    name="nicepy",
    version="0.1",
    author="Gary Chen",
    author_email="gkchen@physics.ucla.edu",
    description="NICE experiment data tools",
    requires=['numpy', 'scipy', 'pint', 'matplotlib', 'pandas']
)
