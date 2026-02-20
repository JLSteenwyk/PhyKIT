import runpy


def test_module_main_invokes_phykit(monkeypatch):
    calls = {"count": 0}

    def fake_phykit():
        calls["count"] += 1

    monkeypatch.setattr("phykit.phykit.Phykit", fake_phykit)

    runpy.run_module("phykit.__main__", run_name="__main__")

    assert calls["count"] == 1
