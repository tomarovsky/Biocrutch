#!/usr/bin/env python3
__author__ = 'tomarovsky'
import cairosvg
import sys

def convert_svg_to_png(input_svg_path, output_png_path, scale_factor=3):
    cairosvg.svg2png(url=input_svg_path, 
                     write_to=output_png_path, 
                     scale=scale_factor)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script_name.py SVG PNG Scale")
        sys.exit(1)

    input_svg_path = sys.argv[1]
    output_png_path = sys.argv[2]
    scale_factor = float(sys.argv[3])

    convert_svg_to_png(input_svg_path, output_png_path, scale_factor)


