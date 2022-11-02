#!/usr/bin/python3
__author__ = 'tomarovsky'
from ete3 import TextFace, Tree, faces, AttrFace, TreeStyle, NodeStyle
from argparse import ArgumentParser


def mylayout(node):
    if node.is_leaf():
        N = AttrFace("name", fgcolor="black", text_prefix="  ", fstyle="italic") # fsize=11
        faces.add_face_to_node(N, node, 0)
        # F = AttrFace("support", fsize=8, fgcolor="green", text_prefix="")
        # faces.add_face_to_node(F, node, 0)
        node.img_style["size"] = 1
        node.img_style["shape"] = "circle"
        node.img_style["fgcolor"] = "Black"
    else:
        node.img_style["shape"] = "circle"
        node.img_style["size"] = 3
        if node.support < 70:
            node.img_style["fgcolor"] = "OrangeRed"
        else:
            node.img_style["fgcolor"] = "LimeGreen"


def main():
    t = Tree(args.input)
    for i in t.get_leaves(): # 'Homo_sapiens' -> 'Homo sapiens'
        i.name = i.name.replace("_", " ")
    if args.outgroup:
        t.set_outgroup(args.outgroup)
    else:
        t.unroot()

    ts = TreeStyle()
    ts.mode = "r"
    ts.layout_fn = mylayout
    ts.show_leaf_name = False

    ts.show_branch_length = True
    ts.show_branch_support = True
    ts.branch_vertical_margin = -12
    t.render(f"{args.output}.length_and_support_tree.svg", w=500, units="px", tree_style=ts)
    ts.show_branch_length = False
    ts.show_branch_support = True
    ts.branch_vertical_margin = -4
    t.render(f"{args.output}.only_support_tree.svg", w=500, units="px", tree_style=ts)
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.branch_vertical_margin = -2
    t.render(f"{args.output}.only_tree.svg", w=500, units="px", tree_style=ts)
    # t.show(tree_style=ts)

if __name__ == "__main__":
    parser = ArgumentParser(description="script to visualize phylogenetic trees using ete3 (required python < 3.10)")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, help="NEWICK file")
    group_required.add_argument('-o', '--output', type=str, help="outfile name")
    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-g', '--outgroup', type=str, default=False, help="outgroup species name (default = unrooted)")
    args = parser.parse_args()
    main()

