#!/usr/bin/python3
__author__ = 'tomarovsky'
from ete3 import TextFace, Tree, faces, AttrFace, TreeStyle, NodeStyle
from argparse import ArgumentParser


def mylayout(node):
    if node.is_leaf():
        N = AttrFace("name", fgcolor="black", text_prefix="  ", fstyle="italic", fsize=12)
        faces.add_face_to_node(N, node, 0)
        # F = AttrFace("support", fsize=8, fgcolor="green", text_prefix="")
        # faces.add_face_to_node(F, node, 0)
        node.img_style["size"] = 1
        node.img_style["shape"] = "circle"
        node.img_style["fgcolor"] = "Black"
        node.dist = 0


def newick_to_nhx(newick_file) -> str:
    with open(newick_file, 'r') as file:
        tree_string = ''
        newick = file.readline().strip().split("'")
        tree_string += newick[0]
        for i in range(1, len(newick), 2):
            line = ''
            flag = True
            for s in newick[i+1]:
                if s == ")" or s == ",":
                    if flag is True:
                        nhx = newick[i].replace(',', '.').replace(';', ':')[1:]
                        line += f"[&&NHX:{nhx}{s}"
                        flag = False
                    else:
                        line += s
                else:
                    line += s
            tree_string += line
        # print(tree_string)
        return tree_string


def main():
    t = Tree(newick_to_nhx(args.input))
    if args.outgroup:
        t.set_outgroup(args.outgroup)
    else:
        t.unroot()
    for i in t.get_leaves(): # 'Homo_sapiens' -> 'Homo sapiens'
        i.name = i.name.replace("_", " ")
    ts = TreeStyle()
    ts.mode = "r"
    ts.layout_fn = mylayout
    ts.show_leaf_name = False

    for n in t.traverse():
        if hasattr(n,"q1"):
            n.add_face(TextFace("  q1={:.2f}  ".format(float(n.q1))), column=2, position="branch-top")
            n.add_face(TextFace("  q2={:.2f}  ".format(float(n.q2))), column=2, position="branch-top")
            n.add_face(TextFace("  pp1={:.2f}  ".format(float(n.pp1))), column=2, position="branch-top")
            n.add_face(TextFace("  pp2={:.2f}  ".format(float(n.pp2))), column=2, position="branch-top")
            n.add_face(TextFace("  EN={:.2f}  ".format(float(n.EN))), column=2, position="branch-top")

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

