# AeroOpt v0.2.0 — Geometria e infraestrutura de solvers

Base geométrica 2D para a futura otimização aerodinâmica multifidelidade de
aerofólios multielementos. A v0.1 importa um elemento principal e um flap
independentes, normaliza, posiciona, valida, exporta e desenha os contornos.
**Não calcula aerodinâmica, CFD, desempenho de voo ou otimização.**

## Interface de solver v0.2.0

O pacote `aeroopt.solvers` introduz condições operacionais imutáveis, contrato
abstrato de adaptadores, identidade científica SHA-256, tentativas isoladas,
resultados tipados, logs e tratamento de falhas. O único adaptador entregue é
`MockSolver`, explicitamente não físico, com `cl`, `cd` e `cm` sempre ausentes.

```powershell
.\.venv\Scripts\python.exe examples/demo_solver_interface.py --output outputs/solver_interface_reproducao
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe -m pytest -q --reference-data=outputs/reference_data_v011
```

A demonstração cobre sete estados, uma repetição da mesma análise e uma rejeição
geométrica: nove tentativas, sem sobrescrita. Os resultados verificados estão
em `outputs/solver_interface_v020_final/`, com JSON, DAT separados e logs.
Nenhum gráfico de coeficientes é produzido. O timeout do mock usa um subprocesso
real bloqueado, encerrado pelo controle de árvore de processos do Windows.

```python
from pathlib import Path
from aeroopt.solvers import MockSolver, OperatingConditions, SolverRunner

runner = SolverRunner(Path.cwd(), "outputs/my_attempts")
conditions = OperatingConditions(alpha_deg=6.0, reynolds=200000.0, mach=0.05)
result = runner.run(MockSolver(), geometry, conditions,
                    min_clearance=0.002, timeout_s=0.5)
```

Nesse trecho, `geometry` é uma montagem válida criada pela API geométrica abaixo.
Reynolds usa a corda principal; alpha está em graus, com direção de velocidade
`(cos(alpha),sin(alpha))` no quadro global. Falha de processo, parsing e
convergência são estados separados. O hash não depende da pasta/horário; cada
tentativa possui UUID próprio. Não há paralelismo de análises nesta versão.

Consulte [SOLVER_INTERFACE.md](docs/SOLVER_INTERFACE.md) para contratos, esquema,
convenções, limites de timeout e implementação de adaptadores reais; veja
[VALIDATION_V020.md](docs/VALIDATION_V020.md) para os testes e a preparação Git.

## Validação geométrica v0.1.1

A v0.1.1 acrescenta geração analítica NACA de quatro dígitos e validação com
NACA 0012, NACA 4412 e S1223 publicados na UIUC. Os originais publicados ficam
em `outputs/`, preservados byte a byte com SHA-256; sua redistribuição não foi
autorizada explicitamente pelas páginas consultadas, portanto não são incluídos
nos arquivos versionáveis. `examples/reference_catalog.json` registra URLs,
hashes revisados, contagens e referências explícitas. Esses arquivos são
coordenadas publicadas, não apresentados como medições experimentais.

Execução offline, usando apenas geometrias **geradas analiticamente**:

```powershell
.\.venv\Scripts\python.exe examples/validate_real_airfoils.py --output outputs/validation_offline
```

Para baixar os três originais fixados por hash e incluí-los no estudo:

```powershell
.\.venv\Scripts\python.exe examples/validate_real_airfoils.py --fetch-uiuc --output outputs/validation_com_uiuc
```

Para reproduzir a execução já realizada, reutilizando os originais locais:

```powershell
.\.venv\Scripts\python.exe -m pytest -q --reference-data=outputs/reference_data_v011
.\.venv\Scripts\python.exe examples/validate_real_airfoils.py --reference-dir outputs/reference_data_v011 --output outputs/validation_reproducao
```

Sempre escolha um `--output` novo/vazio; arquivos anteriores não são sobrescritos.
Se a fonte mudar de hash, a execução exige revisão explícita, preservando o
download recebido. Sem `--reference-data`, somente os 21 casos com arquivos UIUC
são ignorados no pytest; toda a validação analítica permanece offline. Temporários
do pytest ficam automaticamente em uma pasta nova de `.tmp/` quando `--basetemp`
não é fornecido. Nenhuma dependência adicional foi necessária nesta versão.

O exemplo gera perfis normalizados, comparação LE automática/explícita, 32 linhas
de sensibilidade, 45 transformações de normalização quando os três publicados
estão presentes, configurações válidas/rejeitadas, exportações verificadas,
proveniência e tabelas CSV/JSON. A execução final está em `outputs/validation_v011_final/`;
a metodologia e os resultados estão em [VALIDATION_V011.md](docs/VALIDATION_V011.md).

```python
from aeroopt.geometry import Airfoil, naca4_points

points = naca4_points("4412", n_points=161, closed_te=False, spacing="cosine")
foil = Airfoil.from_points(points, leading_edge_index=80)
```

O índice explícito preserva a origem nominal NACA. No perfil cambado, esse ponto
pode diferir do mínimo x e do ponto mais distante da fuga. O modo automático
anterior foi preservado; consulte os erros quantificados no relatório.

As correções de tolerância verificam duplicações e área degenerada na escala
própria de cada elemento. `LENGTH_TOL=1e-10` continua sendo a tolerância global
de contato/folga. O relatório JSON agora informa explicitamente os referenciais
das tolerâncias, mantendo os campos anteriores e o esquema DAT/JSON v0.1.

## Ambiente e instalação

Desenvolvido e verificado com Windows 11 e Python 3.13.15. Execute na raiz do
projeto, utilizando a `.venv` existente:

```powershell
New-Item -ItemType Directory -Force .tmp,.cache\pip | Out-Null
$env:TEMP = Join-Path (Get-Location) '.tmp'
$env:TMP = $env:TEMP
$env:PIP_CACHE_DIR = Join-Path (Get-Location) '.cache\pip'
$env:PIP_DISABLE_PIP_VERSION_CHECK = '1'
.\.venv\Scripts\python.exe -m pip install -e '.[dev]'
.\.venv\Scripts\python.exe -m pytest -q --basetemp=.tmp/pytest
.\.venv\Scripts\python.exe examples/demo_geometry.py
```

Não é necessário ativar ou recriar o ambiente. Não use o Python global. As
dependências são NumPy, Shapely, Matplotlib e pytest (extra `dev`). Os intervalos
de versões estão no `pyproject.toml`; cada exportação registra as versões
efetivas. Para repetir o ambiente validado, instale primeiro
`requirements-dev-lock.txt` e depois o projeto editável. O lock documenta o
ambiente Windows/Python 3.13 verificado, não garante portabilidade entre sistemas.

O exemplo usa backend `Agg`, não abre janelas e salva em `outputs/demo/`:

- `comparison.png` e `comparison.svg`: seis configurações, incluindo rejeições.
- `deflected_detail.png` e `deflected_detail.svg`: detalhe de uma configuração.
- Uma pasta por configuração válida com `main.dat`, `flap.dat`, `metadata.json`.
- `summary.json`: parâmetros, métricas e motivos de rejeição de todos os casos.

O diretório de saída precisa estar vazio para impedir sobrescritas. Para repetir:

```powershell
.\.venv\Scripts\python.exe examples/demo_geometry.py --output outputs/demo_002
```

Os dois perfis em `examples/data/` são polígonos sintéticos de demonstração,
com 9 e 7 pontos; não possuem desempenho aerodinâmico validado.

## Uso da API

```python
from aeroopt.geometry import MultiElementAirfoil

base = MultiElementAirfoil.from_dat(
    main="examples/data/demo_main.dat",
    flap="examples/data/demo_flap.dat",
    flap_chord_ratio=0.30,
    flap_reference_position=(0.95, -0.02),
)
geometry = base.configure(
    flap_deflection_deg=25.0,
    flap_dx=0.08,
    flap_dy=-0.04,
)
report = geometry.validate(min_clearance=0.002)
print(report.to_dict())
if report.valid:
    geometry.export("outputs/candidate_001", min_clearance=0.002)
ax = geometry.plot()
ax.figure.savefig("outputs/candidate_001.png", dpi=180)
```

`configure` **retorna um novo objeto**. Os três valores são absolutos; argumentos
omitidos voltam a zero. O objeto base, os perfis e seus arrays são imutáveis.
`main_points` e `flap_points` estão sempre no referencial global normalizado pela
corda principal. Para mudar razão de cordas ou posição nominal, construa uma
nova base com `from_dat`, ou use `FlapParameters` e o construtor diretamente.

`Airfoil.original_points` preserva os valores importados, sem reordenação;
`normalized_points` fornece o contorno no referencial próprio do perfil.
`leading_edge`, `trailing_edge` e `chord` estão nas unidades originais; `area`
está nessas unidades ao quadrado e `normalized_area` na corda própria ao quadrado.
Áreas de contornos inválidos retornam `None`. `read_dat` apenas lê os pares;
`Airfoil.from_dat` também verifica duplicações e degeneração e normaliza.

Para comparar configurações e salvar as duas extensões:

```python
from aeroopt.geometry.plotting import compare_configurations, save_figure

fig = compare_configurations(
    [base, geometry], labels=["Nominal", "Deslocada"], min_clearance=0.002
)
save_figure(fig, "outputs/comparison")
```

## Formato DAT suportado

Texto UTF-8 (BOM opcional), um contorno ordenado **TE → superfície → LE → outra
superfície → TE**. A primeira linha não vazia pode ser um cabeçalho cujo primeiro
token comece por letra ou `#` e não seja um número (`NACA 0012` é válido). Depois,
cada linha não vazia deve conter exatamente dois números reais separados por
espaços ou tabulações. Notação científica `e`/`E` é aceita. Linhas vazias são
ignoradas. Uma borda de fuga pode ter espessura finita ou repetir o ponto inicial
na última linha. Ambos os sentidos de percurso são aceitos.

```text
Demo profile
1.0 0.01
0.5 1.0e-1
0.0 0.0
0.5 -1.0e-1
1.0 -0.01
```

Não são suportados Lednicer com contagem de pontos, vários blocos/superfícies,
vírgulas, comentários após coordenadas, notação Fortran `D`, colunas extras ou
coordenadas não finitas. Um par inicial de inteiros positivos ambos maiores que
1 é rejeitado como ambíguo com uma tabela de contagens; para coordenadas reais
nessa situação, use `Airfoil.from_points` explicitamente ou declare o primeiro
token do cabeçalho como `AeroOpt_COORDINATES_v0.1`. Esse marcador declara que
todas as linhas seguintes são coordenadas (nunca contagens); é também o
cabeçalho usado nas exportações, permitindo reimportar qualquer posição global.
Erros de parsing
indicam arquivo/linha; erros geométricos indicam o tipo e, para duplicações, os
índices base zero. Não há conversão silenciosa de formatos.

## Validação e exportação

`validate` retorna `ValidationReport`, com `valid`, `issues`, `min_distance`,
`min_clearance`, `overlap_area` e tolerâncias. Cada motivo contém `code`,
`message` e `element`. Contato, cruzamento, contenção sólida, autointerseção,
degeneração, duplicações e folga insuficiente são distinguidos. Métricas entre
elementos inválidos ficam `None`. Entrada escalar inválida gera `ValueError`;
erro de formato gera `AirfoilFormatError`.

`export` revalida com seu próprio argumento explícito `min_clearance` (padrão
zero). Uma chamada anterior a `validate` não altera o objeto e não configura
uma política persistente. Passe a mesma restrição a ambas as funções. Uma
geometria rejeitada gera `InvalidGeometryError`, com o relatório em `.report`,
antes de criar arquivos. Arquivos existentes são protegidos de sobrescrita.

Os DAT exportados contêm 17 algarismos significativos e exatamente os pontos
globais usados na validação, preservando separadamente os elementos. Ao
reimportar, use `read_dat` ou `Airfoil.original_points` para recuperar esse
referencial global. A propriedade `normalized_points` normaliza novamente cada
elemento em seu referencial próprio e, portanto, não é uma reconstrução da
montagem. Use o JSON para recuperar parâmetros e proveniência.

As operações de escrita não constituem uma transação entre os três arquivos:
uma falha de I/O pode deixar uma exportação parcial. O importador preserva os
valores numéricos, cabeçalho e caminho de origem, não uma cópia byte a byte da
formatação textual do arquivo de entrada.

## Documentação e limites

- [Contexto e evolução científica](docs/PROJECT_CONTEXT.md).
- [Convenções, equações, tolerâncias e limites da normalização](docs/GEOMETRY_CONVENTIONS.md).
- [Instruções para próximas sessões](AGENTS.md).
- [Execuções, testes e inspeção dos artefatos](docs/VALIDATION.md).

A borda de ataque é uma referência **amostrada**, não a solução exata da curva
de nariz. A identificação requer o contrato TE-LE-TE; veja o índice explícito
para casos ambíguos. Não há reparo, suavização, interpolação ou redistribuição.
Gap e overlap aerodinâmicos permanecem pendentes de definição: a distância
mínima entre contornos não é apresentada como gap. Toda validação é planar e
poligonal; não inclui tolerâncias de fabricação nem viabilidade aerodinâmica.

