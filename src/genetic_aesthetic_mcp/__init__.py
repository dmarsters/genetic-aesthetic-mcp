"""Genetic Aesthetic Enhancement MCP Server

Provides deterministic codon→amino acid→visual parameter mappings
for image generation prompt enhancement.
"""

__version__ = "0.1.0"

from .server import mcp

__all__ = ["mcp"]
