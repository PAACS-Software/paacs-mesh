# UI bridge examples

Build the optional shared adapter with `-DPAACS_BUILD_UI_BRIDGE=ON`.

```sh
python3 examples/ui/demo_server.py --library /absolute/path/to/libpaacs_ui.so
```

Open the printed URL. The example can generate meshes, inspect entity tags,
cancel work and load/save recipes. Use `recipes/cylinder-tet4.json`,
`cylinder-tet10.json`, `profiled-cylinder.json` or `local-ball.json` to start.

For a frontend that launches a backend process, use:

```sh
python3 examples/ui/paacs_bridge.py --library /absolute/path/to/libpaacs_ui.so
```

The second command reads versioned JSON Lines on stdin and responds on stdout.
Add `--hash` if you need the strict result hash; interactive jobs skip that extra
pass by default. Neither command requires third-party Python packages.
The library path/suffix differs on Windows/macOS and multiconfiguration builds.

See [UI_INTEGRATION.md](../../docs/UI_INTEGRATION.md) for the complete contracts,
control coverage, native C++ options, memory limits and next CAD integration work.
See [VALIDATION.md](../../VALIDATION.md) for tested scope and known gaps.
