from setuptools import setup

setup(
    name="vcf-expression-encoder",
    version="0.0.1",
    entry_points={
        "console_scripts":[
            "vcf-expression-encoder = vcf_expression_encoder:main",
        ]
    },
    install_requires=[
        'vcfpy',
    ],
)
