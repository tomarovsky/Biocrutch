#!/usr/bin/env python3
__author__ = 'tomarovsky'
import cairosvg
import sys
import os


def convert_svg_to_png(input_svg_path, output_png_path, scale_factor=3):
    cairosvg.svg2png(url=input_svg_path, write_to=output_png_path, scale=scale_factor)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python script_name.py SVG PNG Scale')
        sys.exit(1)

    input_svg_path = sys.argv[1]
    output_png_path = os.path.splitext(input_svg_path)[0] + '.png'
    scale_factor = float(sys.argv[2])

    convert_svg_to_png(input_svg_path, output_png_path, scale_factor)
