#!/usr/bin/env python3
__author__ = 'tomarovsky'
from ete3 import TextFace, Tree, faces, AttrFace, TreeStyle, NodeStyle, CircleFace
from argparse import ArgumentParser


def mylayout(node):
    if node.is_leaf():
        N = AttrFace('name', fgcolor='black', text_prefix='  ', fstyle='italic', fsize=12)
        faces.add_face_to_node(N, node, column=0)
        node.img_style['size'] = 1
        node.img_style['shape'] = 'circle'
        node.img_style['fgcolor'] = 'Black'
    else:
        # Создаем кружки поддержки без влияния на длину ветви
        if node.support > 90:
            color = 'LimeGreen'
        elif node.support > 70:
            color = '#008cf0'
        elif node.support > 50:
            color = '#883ac2'
        else:
            color = '#ff0000'

        support_circle = CircleFace(3, color=color, style='circle')
        faces.add_face_to_node(support_circle, node, column=0, position='float')


def add_legend(ts):
    legend_items = [
        (' >90', 'LimeGreen'),
        (' 71-90', '#008cf0'),
        (' 51-70', '#883ac2'),
        (' ≤50', '#ff0000'),
    ]
    for text, color in legend_items:
        circle = CircleFace(5, color=color, style='circle')
        label = TextFace(text, fsize=10, fgcolor='black')
        ts.legend.add_face(circle, column=0)
        ts.legend.add_face(label, column=1)


def main():
    t = Tree(args.input)
    t.sort_descendants()

    for i in t.get_leaves():  # 'Homo_sapiens' -> 'Homo sapiens'
        i.name = i.name.replace('_', ' ').replace('GCA ', 'GCA_').replace('GCF ', 'GCF_').replace("'", '')
    if args.outgroup:
        t.set_outgroup(args.outgroup)

    t.ladderize(direction=1)

    ts = TreeStyle()
    ts.mode = 'r'
    ts.layout_fn = mylayout
    ts.show_leaf_name = False
    for n in t.traverse():
        nstyle = NodeStyle()
        nstyle['fgcolor'] = 'Blue'
        nstyle['size'] = 0
        nstyle['vt_line_width'] = 1
        nstyle['hz_line_width'] = 1
        n.set_style(nstyle)

    add_legend(ts)
    ts.legend_position = (0, 0)

    if args.show:
        t.show(tree_style=ts)

    ts.show_branch_length = True
    ts.show_branch_support = True
    ts.branch_vertical_margin = -12
    t.render(f'{args.output}.length_and_support_tree.svg', w=1500, units='px', tree_style=ts)
    ts.show_branch_length = False
    ts.show_branch_support = True
    ts.branch_vertical_margin = -4
    t.render(f'{args.output}.only_support_tree.svg', w=1500, units='px', tree_style=ts)
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.branch_vertical_margin = 0  # -2
    t.render(f'{args.output}.only_tree.svg', w=3000, units='px', tree_style=ts)
    t.render(f'{args.output}.only_tree.png', w=3000, units='px', tree_style=ts)


if __name__ == '__main__':
    parser = ArgumentParser(description='script to visualize phylogenetic trees using ete3 (required python < 3.10)')
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, help='NEWICK file')
    group_required.add_argument('-o', '--output', type=str, help='outfile name')
    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-g', '--outgroup', type=str, default=False, help='outgroup species name (default = unrooted)')
    group_additional.add_argument('--show', action='store_true', help='option to show tree using GUI')
    args = parser.parse_args()
    main()
