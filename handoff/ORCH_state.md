# Orchestrator state (for resuming)

- Wave 0 done: CONTRACT.md, eigenmap/config.py, RUN_LOG.md; orchestrator drafts of eigenmap/linearize.py, eigenmap/opoints.py, tests/gate1_compare.py (untested, written against the contract).
- Wave 1 agents: A (port) stopped with nothing written -> orchestrator is doing A's job itself (eigenmap/engine_model.py, steady_state.py, tests/test_engine.py, handoff/A_report.md).
  A2 (blind replica, sandbox /home/claude/a2_work, must never see $P/eigenmap) stopped early -> RESUME via SendMessage, then copy /home/claude/a2_work/replica into $P/replica and run `python -m tests.gate1_compare`.
  M (mathcore + anchors) and D0 (runners/docs): check their files/handoff reports; resume if unfinished.
- Next: Gate 1 -> freeze engine -> Wave 2 (L, N, G) -> Run 1 (--fast) + v1.zip -> Run 2 + v2.zip -> auditor -> Wave 3 (E, T) -> Gate 3 -> final.zip.
