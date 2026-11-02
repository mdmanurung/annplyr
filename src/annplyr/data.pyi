from annplyr._expr import Source

obs: Source
var: Source
X: Source
raw: Source

def __getattr__(name: str) -> Source: ...
