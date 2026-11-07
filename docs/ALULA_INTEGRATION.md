# Alula via Node.js — AeroOpt v0.4.0

## Escopo e dependência fixada

O adaptador opcional `AlulaSolver` executa diretamente dois pontos de entrada
do [Alula](https://github.com/kennyjensen/alula), sem navegador, port para Python
ou alteração do Hess–Smith:

| Modo | API JavaScript |
|---|---|
| `inviscid` | `src/inviscid/linear-vortex.js`: `solveInviscid` |
| `viscous_incompressible` | `src/viscous/result.js`: `solveViscousAssembly` |

Commit fixado: **`160f43713a1c6eafb0a4b0f035f5b0bb70b016f8`**, de 2026-10-03.
Foi o HEAD observado durante esta integração. Não se acompanha `master` em
execução. Os modos Euler e Euler/BL acoplado não foram integrados.

`src/aeroopt/solvers/_alula/upstream-lock.json` lista SHA-256 dos **47 arquivos**
selecionados: fechamento transitivo dos imports numéricos, `package.json`,
licenças, avisos e proveniência XFOIL. Esses modos usam somente JavaScript
fixado; não precisam de npm, Playwright, WASM, KLU ou compilação. Não se
instalam UI, analytics, PDFs, presets ou datasets de coordenadas. A rotina
genérica `potential/isentropic.js` é uma dependência interna do BL; chamá-la
com Mach=0 não integra os modos Euler.

O snapshot fica isolado em `.deps/alula/<commit>/`, ignorado pelo Git. O código
upstream permanece byte a byte inalterado; os hashes são conferidos na
construção do adaptador, preparação e antes dos imports no processo Node.
O instalador é explícito e não é invocado pelo solver. Não há download ou
atualização automática numa análise.

## Instalação Windows 11

Execute na raiz do repositório, usando a `.venv` existente:

```powershell
# Com rede, somente fontes fixadas; opcionalmente Node portátil local.
.\.venv\Scripts\python.exe scripts/install_alula.py --download --portable-node

# Alternativa offline: checkout/snapshot que contenha os arquivos originais.
.\.venv\Scripts\python.exe scripts/install_alula.py --source-dir C:\caminho\alula
```

O modo offline verifica cada arquivo contra o lock, sem confiar somente no
nome da pasta ou no SHA declarado pelo checkout. Mudanças locais são rejeitadas;
o instalador não sobrescreve fontes diferentes já instaladas. Uma instalação
interrompida pode ser retomada, mas não passa pela verificação enquanto houver
arquivos ausentes. Downloads originais ficam em `outputs/dependency_sources_v040/`.

O Node portátil exercitado foi **v22.23.4, Windows x64**, obtido da distribuição
oficial, com ZIP e executável fixados no lock:

```text
ZIP SHA-256  cbc74b5e0982c9dbb41525fd8fb5dcd32fbee237863f7829072f4639a1ce5f20
EXE SHA-256  dccc52de28d8df2ee22d2062424f8c563070d95d8ad70c54c4bf01d590dec14f
```

O ZIP foi conferido com `SHASUMS256.txt` oficial, preservado localmente. O
instalador reproduz esses mesmos bytes. Extrai apenas `node.exe`, LICENSE,
README e CHANGELOG em `.deps/node-v22.23.4-win-x64/`; não instala npm ou muda
PATH, registro, configuração global ou a `.venv`. Essa receita portátil é
específica de Windows x64. Para outras plataformas, forneça uma instalação
Node >=22 explicitamente. Outras versões não têm promessa de igualdade bit a
bit; a versão e o hash do executável entram na identidade científica.

O adaptador resolve `node_executable` explícito ou Node no PATH; na ausência
deste, tenta somente a localização portátil fixada. Executa `node --version`
sob `AttemptContext.run_process`, com logs únicos em `.tmp/alula-probes/`.
Dependência ausente ou alterada gera uma exceção descritiva na construção;
ausência/alteração posterior impede a preparação. Não se publica um resultado
físico para representar instalação ausente.

## Uso Python

```python
from pathlib import Path
from aeroopt.solvers import (
    AlulaSolver, AlulaViscousSettings, MomentReference,
    OperatingConditions, SolverRunner,
)

reference = MomentReference(.25, 0, "counterclockwise", "q_inf*c_main^2")
solver = AlulaSolver("inviscid", reference)
runner = SolverRunner(Path.cwd(), "outputs/my_alula_attempts")
result = runner.run(solver, geometry, OperatingConditions(4, 1e6, 0),
                    min_clearance=.002, timeout_s=60)

viscous = AlulaSolver(
    "viscous_incompressible", reference,
    viscous_settings=AlulaViscousSettings(ncrit=9, trips=(1, 1)),
)
```

`geometry` é uma montagem `MultiElementAirfoil` válida com dois elementos.
O contrato existente de `SolverCase` contém principal/flap; **não foi alterado
para aceitar um único elemento**. O bridge standalone aceita 1–6 elementos e
é usado para os testes de um elemento. O adaptador Python aceita a montagem
principal/flap existente, usando seus snapshots globais congelados.

Ncrit, trips e referência de momento são escolhas físicas explícitas. A dupla
de trips é aplicada a todos os elementos, com cada fração baseada na própria
corda do elemento; trips individuais por elemento ainda não são expostos.
Extensões `OperatingConditions.additional` são rejeitadas para evitar ignorar
parâmetros solicitados. Os controles numéricos de `AlulaViscousSettings` são
enviados explicitamente e registrados: 60 iterações Newton, 24 estações por
esteira, comprimento de esteira 2 em cordas principais, tolerâncias 1e-8/1e-6
para equações/esteira, até 40 atualizações de esteira e relaxação .8. Esses
padrões reproduzem as escolhas upstream, não uma qualificação de precisão.
Inicialização é `auto` para multielementos e `native` para um elemento;
Jacobiano `analytic`, sem restart ou warm start de outra análise.

## Geometria, forças e limites de entrada

- Coordenadas são as globais de `SolverCase`, em cordas principais, sem
  regeneração, interpolação, suavização, arredondamento, rotação por alpha ou
  redistribuição de painéis. Principal vem antes do flap. As coordenadas de
  retorno são conferidas contra os pontos fornecidos.
- Aceitam-se apenas contornos **CCW**. Alula pode inverter contornos horários;
  AeroOpt os rejeita para impedir uma reordenação silenciosa. O bridge confere
  também a preparação efetiva do upstream. A validação geométrica/folga comum
  do runner permanece obrigatória.
- Fuga afiada requer endpoints exatamente coincidentes. Endpoints distintos
  declaram `straight_TE_cap`: preservam-se todos os pontos e adiciona-se
  explicitamente o único segmento reto já previsto pela geometria AeroOpt.
  Não se aproximam endpoints para fechar uma fuga quase coincidente. Tampas
  previamente subdivididas com índices de cantos arbitrários não são expostas
  por este adaptador; a declaração é o segmento entre os endpoints fornecidos.
- Fuga finita requer `boundary_condition="streamfunction"` explícito no modo
  invíscido; no viscoso essa formulação é fixa. Isso evita que o upstream
  substitua silenciosamente a condição normal solicitada. A pressão da base é
  modelada pelo Alula a partir de Cp terminal, não é uma medição de recirculação.
  O rótulo `ALULA_FINITE_BASE_RESEARCH` é obrigatório. Não qualifica rankings.
- Exigem-se pelo menos nove pontos fornecidos por elemento, até 600 painéis
  sólidos por elemento e 700 painéis molhados totais. No viscoso são pelo menos
  40 painéis molhados por elemento. A tampa acrescentada conta no limite
  sólido, mas não no orçamento molhado. O bridge aplica os mesmos limites.
- Mach deve ser **exatamente zero**. Inviscid aceita alpha em [-90°,90°];
  viscoso em [-20°,20°], Re em [5e4,5e7], Ncrit em [0,14], trips em (0,1].
  Aceitação de faixa não prova validade física em incidência alta. No invíscido,
  Reynolds permanece somente metadado; no viscoso é baseado em `c_main=1`.
- Alpha tem o mesmo sentido do AeroOpt: `V_inf=(cos(alpha),sin(alpha))`,
  em graus. CL é a projeção perpendicular à corrente, para a esquerda do
  escoamento. CD invíscido é `None`, mesmo se existir `pressureDrag` numérico.
  CD viscoso vem do déficit de momento das esteiras, não da soma indiscriminada
  de diagnósticos de pressão/atrito. CL/CM viscosos são os integrados de pressão
  fornecidos pelo upstream; não se adicionou uma contribuição diferente de atrito.
- O upstream tem **CM positivo horário**. O bridge passa a origem solicitada e
  troca o sinal somente quando o usuário pede `counterclockwise`. CL, CD e CM
  da montagem usam corda principal comum; coeficientes por elemento conservados
  em `upstream.elements` também usam a referência comum, não a corda local do
  flap. Os CM brutos por elemento continuam horários; o CM neutro traduzido tem
  a referência/sentido declarados em `coefficients.moment_reference`.

As convenções foram verificadas por integração independente das pressões,
troca de sentido e transferência para origem de momento não padrão, além da
comparação direta com o upstream. Cp invíscido afiado é por painel; a base
finita integra Cp nodal linear e registra médias por segmento; o modo viscoso
retorna Cp em estações nodais, com dados BL/esteira. Não se interpola tudo para
uma distribuição artificial comum.

## Bridge, estados e elegibilidade

`_alula/bridge.mjs` recebe `--upstream`, `--input` e `--output` como argumentos
separados, sem shell. Os arquivos têm schemas próprios
`aeroopt.alula-input.v1`/`aeroopt.alula-output.v1`. A entrada declara modo,
condições, referência, geometria/fechamento, controles e identidade do runtime.
Campos inesperados, números não finitos, modos incompatíveis e geometria
inadequada são rejeitados. Não há `eval`, código JavaScript da entrada ou
seleção arbitrária do módulo a importar: somente os dois entry points fixados
e dependências conferidas por hash. Isso é integridade de um componente
confiável, não um sandbox para código hostil.

A saída machine-readable é exclusivamente o arquivo JSON, criado sem
sobrescrita. Logs ficam separados em stdout/stderr. Conserva diagnósticos,
Cp, BL quando disponível, avisos, limitações, tempos, amostra de memória,
versões, hashes e coeficientes upstream. Valores diagnósticos não finitos de
uma falha são `null` **com uma lista explícita de caminhos e representações**
(`Infinity`, por exemplo). Não viram zeros. Coeficiente upstream não finito
gera erro; coeficientes provisórios finitos são preservados somente na área
bruta, com status, e não nos coeficientes neutros utilizáveis.

Distinguem-se código de saída do processo, `numerical_convergence`,
`model_admissible`, `coefficient_status` e `eligible`. Para aceitação:

| Modo | Evidência exigida |
|---|---|
| Inviscid | `status=solved`, resíduos linear/Kutta e normal ou streamfunction/interior TE finitos <=1e-9; CL/CM e Cp presentes; CD ausente. |
| Viscoso | `status=solved`, resíduos de equações/esteira dentro das tolerâncias declaradas, espessura relativa máxima <=.1, CL/CD/CM e Cp presentes, sem status provisório. |

O limite .1 é a rejeição ampla de admissibilidade do modelo de camada fina do
upstream; estar abaixo dele **não** demonstra precisão experimental. Um caso
`outside-model` pode ter convergido numericamente, mas permanece inelegível.
Coeficientes utilizáveis ficam todos `None` se falhar qualquer requisito.
O parser Python revalida identidades, hash/imutabilidade da entrada específica,
referência, geometria, disponibilidade, diagnósticos, traduções e limitações;
não aceita somente a flag produzida pelo Node.

No contrato neutro, uma análise inelegível retorna `not_converged` como falha
da política de aceitação, com coeficientes ausentes. O JSON específico preserva
se a causa foi não convergência matemática, modelo fora do domínio, ausência
de diagnósticos ou coeficientes provisórios. `physical_coefficients_available`
do framework indica que há coeficientes de um adaptador real aceito; seu nome
histórico não significa validação física, precisão certificada ou aptidão para
otimização. `aerodynamic_ranking_qualified=false` e `physical_validation=false`
permanecem explícitos. Alula não fornece o certificado `cond_2` do Hess–Smith;
resíduos pequenos do sistema discreto não substituem estudo de condicionamento
ou refinamento. Os gates e identidades HS continuam inalterados.

## Identidade, auditoria e Windows

A identidade inclui commit/lock de fontes, modo, versão/hash do bridge,
versão/hash do binário Node, dependências numéricas, condição de contorno,
referência de momento, controles BL/transição, normalizações e política de
aceitação. Não contém caminhos, timestamps ou UUIDs.
Os três recursos bridge/lock/COPYING têm conversão de fim de linha desativada
em `.gitattributes`, preservando seus bytes/hashes em novos checkouts.
Paths de runtime ficam em `alula_runtime.json`. Cada tentativa copia
bridge/lock/COPYING e preserva
inputs/output/logs sob o inventário do runner. Não se escreve numa tentativa
depois que seu resultado foi retornado.

Todos os processos, inclusive a sondagem de versão, usam
`AttemptContext.run_process`: timeout, logs isolados e encerramento de
descendentes por Job Object no Windows. O processo Node não é persistente.
O Node exige grafia absoluta comum para resolver seu módulo inicial; a
integração remove somente o prefixo operacional Win32 `\\?\` dos argumentos
de caminhos, mantendo os mesmos destinos e os caminhos nativos Python para
I/O. Espaços e Unicode são argumentos separados, sem aspas manuais ou shell.

Foi reproduzido `WinError 267` ao criar processo com cwd muito longo, e erro
Node `EISDIR ... lstat C:` com namespace Win32 no módulo inicial/cwd relativo.
Por isso o adaptador rejeita, na preparação, cwd Windows com mais de **240
unidades UTF-16**, reservando espaço para nomes de artefatos. Use uma raiz
curta do `SolverRunner`. Não se truncam nomes, trocam diretórios silenciosamente
ou alteram contratos do supervisor para contornar isso. Caminhos de rede UNC
não foram exercitados; a conversão de grafia é implementada, sem afirmar suporte
operacional validado a esses destinos.

## Licenças e distribuição

Alula declara **GPL-2.0-or-later**; os dez módulos JavaScript XFOIL/Vibefoil
preservam avisos de Mark Drela (copyright 2000) e a proveniência/adaptações de
Kenny Jensen. As cópias de LICENSE, NOTICE.md, XFOIL README/PROVENANCE/LICENSE
são parte obrigatória do snapshot selecionado. Nenhum arquivo upstream foi
relicenciado, alterado ou colocado no Git AeroOpt. Os wrappers JavaScript
desta integração declaram GPL-2.0-or-later e acompanham `COPYING` integral.

A [licença GPL v2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html) impõe
condições de distribuição, incluindo preservação de avisos/licença e requisitos
de fonte correspondente nas formas aplicáveis. Uso não comercial não cria uma
dispensa dessas condições. Manter um processo Node separado não foi usado como
argumento para remover obrigações de uma distribuição combinada. Antes de
distribuir binários/bundles, revisar o conjunto distribuído e o acesso às fontes
correspondentes; esta etapa publica wrappers, receita e hashes, preservando o
componente opcional em fonte separado. O repositório AeroOpt anterior não tinha
LICENSE geral; esta integração não atribui retroativamente uma licença ao
restante do projeto. A escolha de licença global permanece uma decisão do
titular, e esta documentação não resolve automaticamente seu alcance jurídico.

O aviso upstream enumera KLU, cdt2d e datasets; eles **não estão no fechamento
instalado destes dois modos**. Não se presume que GPL do aplicativo relicencie
coordenadas ou PDFs. Os testes utilizam NACA gerados analiticamente, sem baixar
datasets externos. Node tem sua própria LICENSE, preservada com a distribuição
portátil, que permanece fora do Git.

## Execução e evidência

```powershell
.\.venv\Scripts\python.exe examples/demo_alula.py --output outputs/alula_reproducao
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe -m pytest -q --reference-data=outputs/reference_data_v011 --williams-data=outputs/reference_data_v032
```

O exemplo aceita `--node`/`--alula` absolutos, úteis se lançado de outro cwd,
e `--repeats 1..10`. Preserva casos standalone A/B/D/E e suas chamadas diretas
ao upstream, tentativas de dois elementos, pressão/BL, hashes, snapshots e
custos. Os testes dependentes de Node/Alula são ignorados explicitamente se as
dependências opcionais não existirem; uma instalação com hash alterado falha,
em vez de ser tratada como skip. Resultados, limites e recomendação estão em
[VALIDATION_V040.md](VALIDATION_V040.md).
