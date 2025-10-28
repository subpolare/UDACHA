import os
import argparse
from tqdm import tqdm
from cairosvg import svg2png

def convert_svg_to_png(directory):
    if not os.path.isdir(directory):
        print(f'The provided directory ({director}) does not exist.')
        return

    svg_files = [file for file in os.listdir(directory) if file.endswith('.svg')]

    if not svg_files:
        print(f'No SVG files found in the {directory}.')
        return

    for svg_file in tqdm(svg_files, desc = 'Processing SVG files'):
        svg_path = os.path.join(directory, svg_file)
        png_file = os.path.splitext(svg_file)[0] + '.png'
        png_path = os.path.join(directory, png_file)

        try:
            with open(svg_path, 'rb') as svg_file_content:
                svg2png(file_obj=svg_file_content, write_to=png_path)
        except Exception as e:
            print(f'Error converting {svg_file}: {e}')


parser = argparse.ArgumentParser(description = 'Convert all SVG files in a directory to PNG format.')
parser.add_argument('-d', '--directory', required = True, help = 'Directory containing SVG files.')
args = parser.parse_args() 
convert_svg_to_png(args.directory)
