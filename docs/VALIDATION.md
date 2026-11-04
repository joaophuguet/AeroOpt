# Registro de execução — Geometry Engine v0.1

Execução local em 2026-11-03 (data do projeto, America/Sao_Paulo), Windows 11,
diretório `C:\Users\stitc\Documents\Projetos\AeroOpt`.
Interpretador utilizado em todas as operações Python:
`C:\Users\stitc\Documents\Projetos\AeroOpt\.venv\Scripts\python.exe`.

## Instalação e ambiente efetivo

```powershell
.\.venv\Scripts\python.exe -m pip install -e '.[dev]'
.\.venv\Scripts\python.exe -m pip check
```

A primeira tentativa no sandbox falhou ao buscar `setuptools>=77` no PyPI:
`NewConnectionError`, `WinError 10013` (acesso a soquete proibido).
A repetição com acesso de rede autorizado foi concluída. Dependências foram
instaladas exclusivamente na `.venv`; `TEMP`, `TMP` e `PIP_CACHE_DIR` foram
direcionados a `.tmp/` e `.cache/pip/` dentro do workspace. Não houve recriação
do ambiente nem instalação global.

Versões efetivamente utilizadas: Python 3.13.15, NumPy 2.5.4, Shapely 2.2.0,
Matplotlib 3.11.2, pytest 9.2.0. Dependências transitivas estão fixadas em
`requirements-dev-lock.txt`. O backend de build continua declarado no
`pyproject.toml`; o lock de runtime/testes não fixa o ambiente isolado do build.
`pip check` retornou `No broken requirements found.`

## Testes executados

```powershell
.\.venv\Scripts\python.exe -m pytest -q --basetemp=.tmp/pytest_verified
```

Saída da execução final da suíte:

```text
........................................................................ [ 71%]
.............................                                            [100%]
101 passed in 1.82s
```

**101 executados, 101 aprovados, 0 reprovados, 0 ignorados.** A execução
inicial, antes dos testes adicionais de regressão, teve 93 aprovações e nenhuma
falha. As verificações geométricas usam tolerâncias explícitas e resultados
analíticos; não foram usados valores aerodinâmicos.

| Requisito | Evidência principal |
|---|---|
| Importação válida, cabeçalho e notação científica | `test_airfoil.py::test_valid_dat` |
| Rejeição de formatos inválidos/não finitos | `test_airfoil.py::test_bad_dat_has_explanation` |
| Normalização em escalas e orientações distintas | `test_airfoil.py::test_normalization_is_uniform_rigid` |
| Distâncias rígidas e área proporcional à escala² | `test_transforms.py::test_rigid_distances_and_scaled_area` |
| Rotação positiva horária | `test_transforms.py::test_positive_deflection_is_clockwise` |
| Preservação e imutabilidade dos originais | `test_airfoil.py::test_arrays_are_copied_and_cannot_be_unlocked` |
| Autointerseção e rejeição da exportação | `test_export_plotting.py::test_self_intersection_is_rejected_at_configuration_level` |
| Sobreposição e contenção | `test_validation.py::test_overlap_crossing` e `test_containment_has_nonzero_boundary_distance_but_is_invalid` |
| Distância, folga e contato | `test_validation.py::test_analytic_separation_and_clearance` e testes de tolerância/contato |
| Fuga aberta, área e fechamento | `test_airfoil.py::test_open_te_cap_and_midpoint` e `test_validation.py::test_open_contour_uses_closing_segment` |
| Exportação/reimportação no referencial global | `test_export_plotting.py::test_export_global_frame_roundtrip_and_reproducibility` |
| Reprodutibilidade sem acumular transformações | `test_transforms.py::test_absolute_configurations_are_reproducible` |
| Degeneração e duplicações | `test_airfoil.py::test_degenerate_and_duplicate_inputs` e `test_validation.py::test_invalid_contours_report_structured_reasons` |
| Figuras, legendas, eixos e PNG/SVG | `test_export_plotting.py::test_plot_labels_aspect_and_files` |
| Regressão: coordenadas inteiras globais ambíguas | `test_export_plotting.py::test_global_integer_coordinates_export_reimport` |

## Exemplo executado e saídas

```powershell
.\.venv\Scripts\python.exe examples/demo_geometry.py --output outputs/demo_verified
```

Folga exigida: `0.002 c_main`. Perfis sintéticos de demonstração, com 9 e 7 pontos.

| Candidato | Resultado | Distância mínima / c_main | Motivos de rejeição |
|---|---|---:|---|
| downstream | Válido | 0.10031949 | — |
| deflected | Válido | 0.05660389 | — |
| lowered | Válido | 0.08418208 | — |
| intersecting | Rejeitado | 0.00000000 | boundary_intersection, solid_overlap, insufficient_clearance |
| low_clearance | Rejeitado | 0.00050000 | insufficient_clearance |
| displaced | Válido | 0.07045566 | — |

**4 exportados, 2 rejeitados.** Não há diretório DAT/JSON para os rejeitados.
Os quatro candidatos válidos têm arquivos separados `main.dat`, `flap.dat` e
`metadata.json`. O `summary.json` registra os seis casos, incluindo as rejeições.

Figuras finais em `outputs/demo_verified/`:

- `comparison.png` (2520 × 2052) e `comparison.svg`.
- `deflected_detail.png` (1800 × 720) e `deflected_detail.svg`.

A primeira execução também está preservada em `outputs/demo/`. Os diretórios
de resultados são ignorados pelo Git; o código, dados demonstrativos e documentos
permanecem disponíveis para revisão. Para repetir sem sobrescrever resultados,
escolha outro `--output` vazio.

## Inspeções adicionais

- Os PNGs de comparação e detalhe foram abertos e inspecionados visualmente:
  contornos, referências, legendas e parâmetros estão legíveis, com escala igual.
- Os SVGs finais foram analisados como XML, e os PNGs verificados pelo decoder
  de imagens, sem erros de arquivo.
- Para os quatro candidatos finais, os DAT foram reimportados e comparados
  **exatamente**, ponto a ponto, com os arrays reconstruídos pela API a partir
  dos parâmetros do resumo. Todos coincidiram.
- Cada relatório foi recalculado dos DAT exportados e coincidiu com o JSON.
- `git diff --check` terminou sem erros. A branch continua `main`, com o commit
  inicial `02e9b4c`; nenhum commit, remoto ou histórico foi alterado.

## Limites da evidência

Esses testes verificam contratos geométricos e numéricos em exemplos
representativos, não uma prova para qualquer contorno arbitrário. A identificação
da LE continua amostrada e dependente do contrato TE-LE-TE. Distância não é gap
aerodinâmico. Não houve avaliação de aerodinâmica, validade de solver, capacidade
de carga ou aderência ao regulamento da SAE. Veja as decisões pendentes em
`PROJECT_CONTEXT.md` e as tolerâncias em `GEOMETRY_CONVENTIONS.md`.
