push!(LOAD_PATH, ".")

using Aqua
using Chameleon
Aqua.test_ambiguities([Chameleon])
Aqua.test_all(Chameleon; ambiguities = false, unbound_args = false, deps_compat = false, persistent_tasks = false)
