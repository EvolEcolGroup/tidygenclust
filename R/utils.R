# Functions from other packages, to avoid adding dependencies for single
# calls. fct_inorder_base() is an internal function of tidypopgen (a base R
# version of forcats::fct_inorder()), which creates a factor with the levels
# in their order of first appearance; it is copied from tidypopgen when
# tidygenclust is installed
fct_inorder_base <- utils::getFromNamespace("fct_inorder_base", "tidypopgen")
