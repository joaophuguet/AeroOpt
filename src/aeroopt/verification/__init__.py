"""Opt-in numerical reliability studies, separate from physical solver outputs."""

from .refinement import RefinementRecord, RefinementStudy, RankingPolicy, compare_pair

__all__ = ["RefinementRecord", "RefinementStudy", "RankingPolicy", "compare_pair"]
