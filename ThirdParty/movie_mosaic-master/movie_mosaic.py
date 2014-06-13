import sys
import os
import shutil
import errno
import codecs
import re
import argparse
import collections
import jinja2


# Data for a display table cell.
Cell = collections.namedtuple(
    'Cell',
    'row column rc_address image_filename movie_filename'
    )


def main(argv):

    argparser = argparse.ArgumentParser(
        description="Build an HTML gallery for an array of static images and "
        "associated movie clips.")
    argparser.add_argument('-o', '--output_path', default='./output',
                           help="Location for HTML output "
                           "(default: %(default)s)")
    argparser.add_argument('-c', '--copy', action='store_true',
                           help="Copy images and movies to output folder, "
                           "rather than referencing them in-place "
                           "(default: %(default)s)")
    argparser.add_argument('-m', '--movies-only', action='store_true',
                           help="Force generation of a movie-only mosaic, even "
                           "if images are present in data_path. "
                           "(default: %(default)s)")
    argparser.add_argument('data_path',
                           help="Location of the images and movies")
    args = argparser.parse_args(argv)

    ## Build mappings for rc_address->filename, for images and movies.
    rc_address_pattern = r'r(\d+)c(\d+)'
    image_filename_pattern = re.compile(rc_address_pattern +
                                        r'.*\.(?:png|jpg|jpeg)$')
    movie_filename_pattern = re.compile(rc_address_pattern + r'.*\.mp4$')
    image_filenames = {}
    movie_filenames = {}
    errors = False
    for filename in os.listdir(args.data_path):
        image_match = image_filename_pattern.search(filename)
        movie_match = movie_filename_pattern.search(filename)
        if image_match:
            match = image_match
            filename_map = image_filenames
        elif movie_match:
            match = movie_match
            filename_map = movie_filenames
        else:
            continue
        rc_coords = tuple(map(int, match.groups()))
        if rc_coords in filename_map:
            rc_address = rc_formatter(*rc_coords)
            print >>sys.stderr, ("More than one image/movie file found for "
                                 "location {}:".format(rc_address))
            print >>sys.stderr, " ", filename_map[rc_coords], filename
            errors = True
        else:
            filename_map[rc_coords] = filename
    if errors:
        # Exit early if any duplicate files were found.
        return 1

    ## Build a list of Cell objects.
    cells = []
    mkdirp(args.output_path)
    if args.movies_only:
        filenames = movie_filenames
    else:
        filenames = image_filenames
    for (row, column), filename in filenames.items():
        rc_address = rc_formatter(row, column)
        if args.movies_only:
            image_filename = None
            movie_filename = filename
        else:
            image_filename = filename
            try:
                movie_filename = movie_filenames[(row, column)]
            except KeyError:
                movie_filename = None
                print >>sys.stderr, "WARNING: Movie file missing for location", \
                    rc_address, "\n"
        cells.append(Cell(row, column, rc_address, image_filename, movie_filename))

    if args.copy:
        ## Copy images and movies to the output path.
        for cell in cells:
            for filename in cell.image_filename, cell.movie_filename:
                if filename is None:
                    continue
                src_path = os.path.join(args.data_path, filename)
                dest_path = os.path.join(args.output_path, filename)
                shutil.copyfile(src_path, dest_path)
    else:
        ## Transform filenames into absolute URLs.
        cells = [Cell(cell.row, cell.column, cell.rc_address,
                      absurl(cell.image_filename, args.data_path),
                      absurl(cell.movie_filename, args.data_path))
                 for cell in cells]

    ## Build the display table as a dict-of-dicts of Cells.
    row_min = min(c.row for c in cells)
    row_max = max(c.row for c in cells)
    column_min = min(c.column for c in cells)
    column_max = max(c.column for c in cells)
    table = collections.defaultdict(dict)
    for cell in cells:
        table[cell.row][cell.column] = cell
    row_numbers = range(row_min, row_max + 1)
    column_numbers = range(column_min, column_max + 1)

    ## Render the template to a file.
    src_path = os.path.dirname(__file__)
    template_path = os.path.join(src_path, 'templates')
    template_loader = jinja2.FileSystemLoader(template_path)
    template_env = jinja2.Environment(loader=template_loader,
                                      undefined=jinja2.StrictUndefined)
    template = template_env.get_template('index.html')
    data_names = 'table cells row_numbers column_numbers'.split()
    # Take a "slice" of locals corresponding to the names in data_names.
    data = dict(zip(data_names, map(locals().get, data_names)))
    data['movies_only'] = args.movies_only
    content = template.render(data)
    html_filename = os.path.join(args.output_path, 'index.html')
    with codecs.open(html_filename, 'w', 'utf-8') as out_file:
        out_file.write(content)

    ## Copy the static resources to the output path.
    static_output_path = os.path.join(args.output_path, 'static')
    shutil.rmtree(static_output_path, ignore_errors=True)
    shutil.copytree(os.path.join(src_path, 'static'), static_output_path)

    ## Calculate the total size of the output directory.
    output_size = 0
    for root, dirs, files in os.walk(args.output_path):
        for name in files:
            output_size += os.path.getsize(os.path.join(root, name))

    ## Report some useful information.
    print "Data path: {}".format(os.path.abspath(args.data_path))
    print "  rows observed: {}-{}".format(row_min, row_max)
    print "  columns observed: {}-{}".format(column_min, column_max)
    print "  total locations: {}".format(len(cells))
    print
    print "Output path: {}".format(os.path.abspath(args.output_path))
    print "  size: {:.1f} MB".format(output_size/1e6)
    print
    print "Media copied: {}".format('yes' if args.copy else 'no')
    print "Mode: {}".format('movies only' if args.movies_only
                            else 'images and movies')


def rc_formatter(row, column):
    return 'r{}c{}'.format(int(row), int(column))


def mkdirp(path):
    try:
        os.makedirs(path)
    except OSError, e:
        if e.errno != errno.EEXIST: raise


def absurl(path, base):
    if path is not None:
        return 'file://' + os.path.abspath(os.path.join(base, path))
    else:
        return None

    
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
