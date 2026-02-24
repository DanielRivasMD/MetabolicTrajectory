####################################################################################################

module Paths

####################################################################################################

using TOML

####################################################################################################

const CONFIG = TOML.parsefile("src/config/paths.toml")
const STRUCT = CONFIG["data"]
const ROOTS = Dict{Symbol,String}()

####################################################################################################

function register(name::Symbol, path::String)
  ROOTS[name] = path
end

function resolve(node, parent_path)
  children = Dict{String,Any}()

  for child in get(node, "children", [])
    child_node = get(node, child, Dict())
    child_path = joinpath(parent_path, child)
    children[child] = resolve(child_node, child_path)
  end

  return (path = parent_path, children = children)
end

function build(name::Symbol)
  root = ROOTS[name]
  return resolve(STRUCT, root)
end

function path(name::Symbol, key::String)
  tree = build(name)
  parts = split(key, ".")
  node = tree
  for p in parts
    node = node.children[p]
  end
  return node.path
end

function ensure(root)
  mkpath(root.path)
  for (_, child) in root.children
    ensure(child)
  end
end

function mirror_tree(src, dst)
  mkpath(dst.path)
  for (name, src_child) in src.children
    dst_child_path = joinpath(dst.path, name)
    mkpath(dst_child_path)
    dst_child = (path = dst_child_path, children = Dict{String,Any}())
    mirror_tree(src_child, dst_child)
  end
end

function mirror(pair::Pair{Symbol,Symbol})
  src_name, dst_name = pair
  src_tree = build(src_name)
  dst_tree = build(dst_name)
  mirror_tree(src_tree, dst_tree)
end

####################################################################################################

end

####################################################################################################
