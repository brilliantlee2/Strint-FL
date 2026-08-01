import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


class ReleaseLayoutTests(unittest.TestCase):
    def test_release_manifests_exist(self):
        for relative_path in [
            ".gitignore",
            "Cargo.toml",
            "Cargo.lock",
            "environment.yml",
            "requirements.txt",
            "README.md",
            "README_zh-CN.md",
            "run_all.sh",
            "run_all_mixed_species.sh",
            "vendor/rust-htslib/Cargo.toml",
        ]:
            self.assertTrue((PROJECT_ROOT / relative_path).is_file(), relative_path)

    def test_readmes_have_language_switches_and_core_commands(self):
        english = (PROJECT_ROOT / "README.md").read_text(encoding="utf-8")
        chinese = (PROJECT_ROOT / "README_zh-CN.md").read_text(encoding="utf-8")

        for readme in [english, chinese]:
            self.assertIn("README.md", readme)
            self.assertIn("README_zh-CN.md", readme)
            self.assertIn("conda env create -f environment.yml", readme)
            self.assertIn("cargo build --release", readme)
            self.assertIn("bash run_all.sh", readme)
            self.assertIn("bash run_all_mixed_species.sh", readme)

    def test_environment_uses_python_311_and_scanpy_dependencies(self):
        environment = (PROJECT_ROOT / "environment.yml").read_text(encoding="utf-8")
        requirements = (PROJECT_ROOT / "requirements.txt").read_text(encoding="utf-8")
        english = (PROJECT_ROOT / "README.md").read_text(encoding="utf-8")
        chinese = (PROJECT_ROOT / "README_zh-CN.md").read_text(encoding="utf-8")

        self.assertIn("python=3.11", environment)
        self.assertNotIn("python-gil", environment)
        self.assertIn("- nodefaults", environment)
        for dependency in [
            "scanpy>=1.11,<1.12",
            "numba>=0.66,<0.67",
            "python-igraph",
            "leidenalg",
        ]:
            self.assertIn(f"- {dependency}", environment)
        for dependency in ["scanpy>=1.11,<1.12", "numba>=0.66,<0.67", "igraph", "leidenalg"]:
            self.assertIn(dependency, requirements)
        for readme in [english, chinese]:
            self.assertIn("Python 3.11", readme)
            self.assertNotIn("3.14t", readme)
        self.assertRegex(environment, r"(?m)^\s*- rust(?:[=<>]|\s*$)")
        self.assertNotRegex(environment, r"(?m)^\s*- cargo(?:[=<>]|\s*$)")
        for dependency in ["samtools", "minimap2", "bedtools", "pip"]:
            self.assertIn(f"- {dependency}", environment)

    def test_large_runtime_data_are_ignored(self):
        gitignore = (PROJECT_ROOT / ".gitignore").read_text(encoding="utf-8")

        for pattern in ["target/", "*.fastq.gz", "*.bam", "report_new/", "report_new_2/", "vendor.zip"]:
            self.assertIn(pattern, gitignore)

    def test_release_version_is_consistently_0035(self):
        cargo_toml = (PROJECT_ROOT / "Cargo.toml").read_text(encoding="utf-8")
        cargo_lock = (PROJECT_ROOT / "Cargo.lock").read_text(encoding="utf-8")
        runners = "\n".join(
            (PROJECT_ROOT / name).read_text(encoding="utf-8")
            for name in ["run_all.sh", "run_all_mixed_species.sh"]
        )

        self.assertRegex(cargo_toml, r'(?m)^version = "0\.0\.3\+5"$')
        self.assertIn('name = "strint-rust"\nversion = "0.0.3+5"', cargo_lock)
        self.assertIn("StrintRust0.0.3.5", runners)
        self.assertNotIn("StrintRust0.0.2", runners)


if __name__ == "__main__":
    unittest.main()
