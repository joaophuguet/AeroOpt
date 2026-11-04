# Validação geométrica v0.1.1 — 2026-11-04

## Escopo, ambiente e execução

Validação local em Windows 11 / Python 3.13.15 com o interpretador obrigatório
`C:\Users\stitc\Documents\Projetos\AeroOpt\.venv\Scripts\python.exe`.
Pacote atualizado e instalado em modo editável como **0.1.1**, sem novas
dependências. NumPy 2.5.4, Shapely 2.2.0 / GEOS 3.14.1, Matplotlib 3.11.2 e
pytest 9.2.0; `pip check` não encontrou incompatibilidades. A instalação editável
foi atualizada com `pip install -e '.[dev]' --no-deps`, usando cache/temporários
no workspace. Resultados anteriores foram preservados.

Esta é uma validação exclusivamente geométrica. Não se executou solver,
painelamento aerodinâmico, CFD, MSES, RANS, XFOIL, NSGA-II ou estimativa de
coeficientes. Contornos não interceptantes não demonstram eficácia de uma
configuração single-slotted flap nem conformidade com requisitos de voo.

```powershell
.\.venv\Scripts\python.exe -m pytest -q --reference-data=outputs/reference_data_v011 --basetemp=.tmp/v011_reviewed
.\.venv\Scripts\python.exe examples/validate_real_airfoils.py --reference-dir outputs/reference_data_v011 --output outputs/validation_v011_final
```

Resultados finais: **197 testes aprovados, 0 reprovados, 0 ignorados, em 8.77 s**.
São os 101 anteriores e 96 casos novos. Antes das alterações, os 101 testes da
v0.1 passaram em 2.71 s. A execução inicial da nova suíte, antes dos testes do
fluxo offline/download, aprovou 193 casos. Os quatro casos adicionados depois
também passaram na execução final.

Para repetir, escolha outro `--output` vazio e outro `--basetemp`, ou omita
`--basetemp` para que os testes criem uma pasta nova dentro de `.tmp/`.
Os dados UIUC não são necessários aos testes analíticos: sem `--reference-data`,
os 21 casos que dependem dos originais são explicitamente ignorados.
A execução adicional de `.\.venv\Scripts\python.exe -m pytest -q`, sem opções,
confirmou **176 aprovados e 21 ignorados em 8.17 s**, sem acesso à rede e com
temporários criados automaticamente dentro de `.tmp/`.

## Conjuntos de referência e preservação

Os arquivos foram obtidos em **2026-11-04**, preservados sem conversão de
quebras de linha, números ou cabeçalhos. A fonte é a
[UIUC Airfoil Coordinates Database](https://m-selig.ae.illinois.edu/ads/coord_database.html),
no diretório explicitamente ordenado em formato Selig: fuga superior → nariz →
fuga inferior. As unidades são nominais `x/c,y/c`; não são unidades físicas.
Não se afirma que os arquivos sejam medições experimentais de coordenadas.

| Perfil | Pontos | Bytes | Fonte original | Referência explícita (índice base zero) |
|---|---:|---:|---|---|
| NACA 0012 | 131 | 2772 | [n0012.dat](https://m-selig.ae.illinois.edu/ads/coord_seligFmt/n0012.dat) | 65: (0,0) |
| NACA 4412 | 35 | 693 | [naca4412.dat](https://m-selig.ae.illinois.edu/ads/coord_seligFmt/naca4412.dat) | 17: (0,0) |
| S1223 | 81 | 1872 | [s1223.dat](https://m-selig.ae.illinois.edu/ads/coord_seligFmt/s1223.dat) | 45: (0.00005,0.00178) |

Hashes SHA-256 dos bytes originais:

```text
n0012.dat   ac83dca5d32bc4e4db2131820e87981bb71903cfeebee8cefd0818dedbe2b012
naca4412.dat 5c86886f667266f29a5f3a0a48505b9b25ea10539d8a112486c62bf3a5bbd6ad
s1223.dat   f2e06d165d6244d2fc1c7dde9e46c46c4d42a70b183620ea02fd693fbfa5cf53
```

O S1223 foi escolhido como terceiro perfil porque a literatura dos próprios
autores o documenta para baixa velocidade/baixo Reynolds; veja
[Selig e Guglielmo, 1997](https://m-selig.ae.illinois.edu/pubs/GuglielmoSelig-1997-JofAC-S1223.pdf).
Essa contextualização não foi usada como resultado aerodinâmico do AeroOpt.

Foram consultadas as páginas da base, a
[FAQ](https://m-selig.ae.illinois.edu/ads_faq.html) e as
[condições dos dados de túnel](https://m-selig.ae.illinois.edu/pd.html).
As condições dessa última página tratam de dados de desempenho; não se presumiu
que licenciem todos os arquivos de coordenadas. Como não foi estabelecida uma
autorização explícita de redistribuição destes três arquivos, seus bytes e
derivados permanecem em `outputs/`, ignorado pelo Git. URLs, hashes, contagens e
seleções de referência estão em `examples/reference_catalog.json`. Não houve
contato com autores nem envio de mensagens. Uma revisão de hash upstream exige
ação explícita, nunca atualização silenciosa das referências.

Os originais estão em `outputs/reference_data_v011/` e suas cópias byte a byte
em `outputs/validation_v011_final/originals/`, acompanhadas de `provenance.json`.
O catálogo fixa a revisão usada; novas obtenções registram a data efetiva.

## Geração analítica complementar e referências independentes

O gerador `naca4_points` e a avaliação `naca4_surfaces` implementam as equações
NACA de quatro dígitos descritas no
[PDAS — espessura](https://www.pdas.com/naca456thick4.html) e
[PDAS — linha média](https://www.pdas.com/naca456mean2.html).
As equações e as convenções de referência estão em `GEOMETRY_CONVENTIONS.md`.
Todas as geometrias produzidas são identificadas como **geradas analiticamente**,
sem origem experimental. A fuga clássica finita usa `-0.1015`; a alternativa
fechada modifica explicitamente esse coeficiente para `-0.1036`.

O estudo usa N=41,81,161,321, amostragem uniforme e cossenoidal, nos dois perfis
e nos dois tipos de fuga: **32 linhas de sensibilidade**. O gerador não toca
nos arquivos importados. Seus DAT, parâmetros e hashes estão em `analytic/`
e no manifesto de proveniência do resultado.

A correção não é verificada apenas chamando o mesmo gerador outra vez:

- Uma ordenada simétrica é calculada com Decimal em precisão de 50 dígitos.
- Ordenadas publicadas do NACA 4412 em x=0.4 são verificadas com tolerância
  compatível com a tabela de quatro casas decimais (`5e-5 c`).
- Para NACA 0012, a área contínua vem da integral exata do polinômio de espessura:
  `A = 0.082210 c²` (fuga finita) ou `0.081706 c²` (modificada fechada).
- `examples/_naca_reference.py` não importa nem chama o gerador. Usa equações
  escalares expandidas, derivadas analíticas e quadratura Gauss-Legendre após
  `x=u²`. O ponto de troca de linha média é integrado separadamente.

Para justificar a integração: com `h` semiespessura e `g=dyc/dx`, o elemento de
área da faixa normal simétrica é `2h sqrt(1+g²) dx`. O termo ímpar de curvatura
cancela entre as duas faces. As velocidades paramétricas em x têm módulos
`sqrt((s ∓ h*g'/s²)² + h_x²)`, com `s=sqrt(1+g²)`; o perímetro inclui a tampa
de comprimento `2h(1)`. A troca para u remove a singularidade da derivada de
sqrt(x) no nariz. Essas referências se limitam às famílias ensaiadas, sem
reivindicar validade para faixas normais que se dobram ou autointersectam.

Referências contínuas utilizadas:

| Perfil / fuga | Área / c² | Perímetro fechado / c |
|---|---:|---:|
| 0012 finita | 0.0822100000000000 | 2.041791633043174 |
| 0012 fechada | 0.0817060000000000 | 2.039548931156012 |
| 4412 finita | 0.0824974369712582 | 2.050286414458205 |
| 4412 fechada | 0.0819908574995590 | 2.048042311470386 |

A diferença entre ordens 128 e 256 de quadratura foi <=`5.0e-16 c²` para área
e <=`1.34e-15 c` para perímetro; a área simétrica integrada concordou com a
expressão exata em `2.78e-17 c²`. A diferença entre ordens é um diagnóstico,
não um limite rigoroso de erro ou uma prova de exatidão do perímetro.

## Normalização e referências LE/TE

Os cinco conjuntos de base (dois analíticos de 161 pontos e três publicados)
foram transladados, girados em -75°, 0° e 138° e escalados em 1e-6, 1 e 1e6:
**45 casos de normalização**. Os testes também exercitam escalas 1e-9 e 1e9,
usando translação proporcional à escala para não destruir a precisão de entrada.
Referências explícitas foram mantidas durante as mudanças de quadro.

O estudo impõe tolerância de comparação `2e-12`, distinta da política de
contato. Os maiores erros observados foram:

| Métrica | Maior erro |
|---|---:|
| Posição normalizada de uma amostra | 9.2551e-16 c |
| Distância entre amostras, após remover escala | 1.2213e-15 c |
| Área após remover escala², erro relativo | 1.3323e-15 |

A orientação física foi preservada em todos os casos. Arrays de entrada foram
comparados exatamente e os hashes dos originais publicados permaneceram iguais.

Nos três publicados, o índice automático coincidiu com o explicitamente
selecionado. Isso é uma observação destes conjuntos, não garantia geral. Áreas
normalizadas dos publicados: NACA 0012=`0.0821784911651800`,
NACA 4412=`0.0821112500000000`, S1223=`0.0649145848208250`.
Não se atribuem essas diferenças em relação a famílias analíticas apenas a erro
de implementação: discretização, arredondamento e convenção das tabelas diferem.

O NACA 4412 analítico é um contraexemplo à equivalência automática/nominal:
N=161/321 cossenoidal resulta em diferença de referência `0.00349061496976533 c`
e diferença relativa de corda `0.000300138689133478`. Em N=41/81 a amostra mais
distante ainda é a origem nominal. Isso muda a normalização por convenção,
embora cada transformação continue preservando a forma. O índice nominal
explícito evita introduzir essa mudança no estudo de discretização.

Para S1223, a origem nominal exata não está amostrada: a referência selecionada
é `(0.00005,0.00178)`, e a corda resultante é `0.9999515842779589` nas unidades
publicadas. Não se declara erro zero relativo à LE contínua desconhecida.

## Resultados quantitativos de discretização

Tabela abaixo: fuga finita, estações cossenoidais, erros absolutos contra as
referências contínuas independentes (não contra a malha mais fina).

| Perfil | N | Erro de área / c² | Erro de perímetro / c | Erro de distância / c |
|---|---:|---:|---:|---:|
| 0012 | 41 | 3.32400e-4 | 1.15163e-3 | 8.67e-18 |
| 0012 | 81 | 8.31926e-5 | 2.87294e-4 | 8.67e-18 |
| 0012 | 161 | 2.08039e-5 | 7.17577e-5 | 8.67e-18 |
| 0012 | 321 | 5.20134e-6 | 1.79357e-5 | 8.67e-18 |
| 4412 | 41 | 3.42154e-4 | 1.21094e-3 | 9.00744e-5 |
| 4412 | 81 | 8.56404e-5 | 3.02072e-4 | 9.00744e-5 |
| 4412 | 161 | 2.14206e-5 | 7.54750e-5 | 1.84691e-6 |
| 4412 | 321 | 5.35553e-6 | 1.88659e-5 | 1.84691e-6 |

Com N=321 uniforme, os erros de área/perímetro foram, respectivamente,
`3.70176e-5 / 8.00979e-4` para 0012 e `3.77384e-5 / 8.23197e-4` para 4412.
A distribuição cossenoidal resolve o nariz melhor neste estudo. Não se recomenda
uma redistribuição automática de arquivos existentes com base nisso.

Para distância, o principal é sempre NACA 0012 de fuga finita, e o flap tem
razão 0.3, posição nominal `(1.01,0)`, deflexão zero e referência nominal
explícita. O limite a jusante do sólido principal é x=1. O mínimo x contínuo
do flap projeta-se dentro da tampa vertical da fuga principal; a distância
global é, portanto, `0.01 + 0.3 * min(x_flap_local)`.

Para 0012, esse valor é exatamente `0.01 c`. Para 4412, resolve-se a derivada
analítica `dx_upper/du=0` por bisseção em `[0,0.25]`; o mínimo resultante tem
`x=-0.000300247899390`, `y=0.003037524455609`. O valor global da distância é
`0.00990992563018304 c`, com projeção `0.3*y=0.0009112573` dentro da tampa
de semiespessura 0.00126. O resíduo da derivada ficou abaixo de `6e-18`.
Essa construção é específica desses perfis/posição, não um algoritmo geral
de distância entre curvas contínuas arbitrárias.

**N=161 e N=321 produziram a mesma distância poligonal 0.00991177253970510 c
no caso 4412**, com erro não nulo `1.84691e-6 c`. As duas malhas contêm a mesma
amostra extrema relevante. Logo, concordância entre discretizações não prova
exatidão. Com estações uniformes, esse nariz não é resolvido nem em N=321:
distância 0.01 c, erro `9.00744e-5 c`. Os 32 casos de sensibilidade foram
geometricamente válidos para folga requerida 0.002 c; isso não torna seus erros
de discretização intercambiáveis.

## Defeitos reproduzidos e correções

1. **Falsa degeneração por escala:** na v0.1, um flap 0012 de razão 1e-6 ou
   1e-7 era rejeitado por `Closed area <= area tolerance`, apesar de a forma
   relativa ser a mesma e os pontos permanecerem resolvidos em float64.
2. **Falsa duplicação por discretização:** em um flap de razão 1e-5, inserir
   um ponto ao longo de um segmento já existente mudava o status para
   `duplicate_points`, mesmo sem alterar a região sólida.

Antes da correção, os quatro testes de `test_tolerance_regressions.py` deram
**3 falhas e 1 aprovação**. Depois, os testes de regressão e os testes antigos
de validação deram **21 aprovações**. A suíte final manteve os 101 anteriores.

As constantes `LENGTH_TOL=1e-10` e `AREA_TOL=1e-12` foram preservadas. A correção
usa distâncias/D e área/D² para qualidade intrínseca, onde D é a diagonal da
caixa do próprio elemento. Essas razões são invariantes a escala uniforme.
Contato/folga continuam absolutos em `c_main`. O relatório inclui agora
`element_tolerance_frame`, `area_tolerance_units`, `duplicate_tolerance_units`
e `contact_clearance_tolerance_units`. Campos anteriores, API de configuração,
marcador DAT e identificador de esquema JSON foram preservados. Há uma mudança
intencional da interpretação da tolerância de área no relatório; consumidores
devem ler os novos campos de unidade.

A referência automática diferente da nominal foi registrada como **limitação
comprovada do método escolhido**, não corrigida por uma mudança silenciosa de
LE. O estudo usa a opção de índice explícito já disponível na v0.1.

## Configurações difíceis e round trips

Foram avaliados 10 casos analíticos e 3 combinações dos perfis publicados:
**13 configurações, 8 válidas e 5 rejeitadas**. Deflexões 0°, +25° e -20°,
razões 0.2/0.3/0.03/1e-7, fuga finita, contato, cruzamento e contenção foram
exercitados. Os testes publicados também incluem razão 0.4.

| Caso | Distância / c | Resultado |
|---|---:|---|
| zero | 0.0200000000 | Válido |
| positive | 0.0427520880 | Válido |
| negative | 0.0426880504 | Válido |
| insufficient | 0.0000500000 | Folga insuficiente para 0.002 c |
| contact | 0 | Contato de contornos, área de sobreposição zero |
| crossing | 0 | Cruzamento e sobreposição sólida |
| contained | 0.0533267324 | Sobreposição sólida, sem cruzamento de contornos |
| above_tolerance | 2.00000017e-10 | Válido, acima do limite 1e-10 |
| below_tolerance | 5.00000041e-11 | near_contact |
| small_resolved_flap | 0.0200000000 | Válido após correção intrínseca |
| published_1 | 0.1072417251 | 0012 principal / 4412 flap, válido |
| published_2 | 0.1072188006 | 4412 principal / S1223 flap, válido |
| published_3 | 0.1064463564 | S1223 principal / 0012 flap, válido |

Nas oito exportações válidas, a reimportação de `main.dat` e `flap.dat` deu
**erro máximo exatamente zero** nos valores float64. Distância e relatório
recalculados coincidiram com o JSON; identificação, nomes de elementos e
parâmetros foram preservados. Há fuga de espessura finita nos round trips.
A reimportação compara o referencial global por `read_dat`, não uma nova
normalização local do flap. Resultados rejeitados têm relatório, sem exportação
de contornos como se fossem válidos.

## Artefatos e testes novos

`outputs/validation_v011_final/` contém:

- `normalized_profiles.png/.svg`: LE automática e explícita nos cinco perfis.
- `discretization_errors.png/.svg`: erros de área, perímetro e distância.
- `configurations_1.png/.svg` até `configurations_3.png/.svg`: configurações e
  rejeições, com a folga efetiva individual indicada no título.
- `reference_metrics.csv/.json`, `normalization_invariance.csv/.json`,
  `discretization.csv/.json`, `roundtrips.csv/.json`.
- `continuous_references.json`, `configuration_results.json`, `summary.json`,
  `provenance.json`, `analytic/`, `originals/`, `exports/`.

As figuras foram inspecionadas visualmente, e os arquivos PNG/SVG finais
verificados pelo decoder de imagens/parser XML. As versões intermediárias
permanecem em `outputs/validation_v011_first/` (antes da atualização de versão)
e `outputs/validation_v011/` (antes de vincular os DAT analíticos à proveniência
de cada exportação). Os resultados finais usam 0.1.1 e todos os elementos
exportados possuem arquivo de origem registrado.

Testes novos: `test_naca.py`, `test_reference_validation.py`,
`test_published_airfoils.py`, `test_tolerance_regressions.py` e
`test_validation_workflow.py`. Cobrem fórmulas, referências independentes,
semelhanças, preservação, stress multielemento, discretização, round trips,
hashes, download simulado, execução offline e proteção contra sobrescrita.
Testes de downloads não acessam rede: usam bytes sintéticos declarados.

## Prontidão para integração e limites remanescentes

O núcleo está pronto para **iniciar um adaptador de solver em escopo controlado**:
montagens são determinísticas, contornos separados, dados preservados, erros
estruturados e round trips exatos nos casos exercitados. A compatibilidade
geométrica da v0.1 foi mantida, com as duas correções justificadas.

Não há evidência suficiente para confiar diretamente em uma campanha de
otimização. Antes disso, é necessário revisar LE/TE de cada entrada, definir
resolução e qualidade de painelamento por solver, validar suas convenções de
fuga, formalizar gap/overlap e exercitar condições de falha/convergência. A
distância mínima global permanece distinta de gap aerodinâmico. Exportação
ainda não é transacional entre os três arquivos, e a reconstrução de montagem
por JSON não está implementada.

Escalas ou translações extremas podem perder informação em float64; as
tolerâncias intrínsecas não recuperam pontos já colapsados. As referências
contínuas não cobrem todo NACA de quatro dígitos ou todo perfil da UIUC. Não há
certificado rigoroso de erro para quadratura/perímetro, nem para a LE contínua
do S1223. As condições de redistribuição das coordenadas externas continuam
pendentes, e sua disponibilidade futura deve ser tratada pelo catálogo de hashes.
