import Symbolics: unwrap, wrap
import ModelingToolkit: VariableConnectType

# =============================================================================
# BG types

struct bg end
struct bgeffort end
struct bgflow end
struct op end
struct tpgy end
struct tptf end
struct j0 end
struct j1 end

# BG type functions

set_bg_metadata(s, type) = wrap(setmetadata(unwrap(s), bg, type))
get_bg_junction(s) = getmetadata(unwrap(s), bg)
update_mtk_con(s, type) = wrap(setmetadata(unwrap(s), VariableConnectType, type))
