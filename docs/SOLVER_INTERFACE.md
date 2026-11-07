# Solver Interface — AeroOpt v0.2.0

> Contratos preservados na v0.3.0. O escopo abaixo descreve a entrega histórica
> v0.2.0; o adaptador físico adicionado posteriormente é documentado em
> [PANEL_METHOD.md](PANEL_METHOD.md). Não houve mudança dos esquemas de entrada,
> resultado ou das quatro etapas do adaptador. Cp/circulação/resíduos são
> artefatos específicos inventariados, sem campos obrigatórios para futuros
> adaptadores MSES/SU2.

> A v0.4.0 preserva esses mesmos contratos e acrescenta o adaptador externo
> [Alula via Node.js](ALULA_INTEGRATION.md). Os estados detalhados de convergência,
> admissibilidade e coeficientes provisórios ficam no artefato específico;
> análises inelegíveis mantêm coeficientes neutros ausentes. O runner continua
> limitado à montagem principal/flap, enquanto o bridge standalone aceita
> casos de um elemento para verificação. Não houve alteração dos schemas v1.

## Contrato e escopo

Esta versão fornece contratos e execução auditável de análises 2D. O único
adaptador entregue é `MockSolver`, explicitamente **não físico**. Ele não prevê
forças, momentos, polares, estol nem desempenho. A interface permite adaptadores
futuros por processo externo ou biblioteca Python, sem depender de formatos
MSES/XFOIL/RANS. O pacote geométrico e as correções da v0.1.1 são preservados.

```python
from pathlib import Path
from aeroopt.geometry import MultiElementAirfoil
from aeroopt.solvers import MockSolver, OperatingConditions, SolverRunner, SolverResult

geometry = MultiElementAirfoil.from_dat(
    "examples/data/demo_main.dat", "examples/data/demo_flap.dat",
    flap_chord_ratio=0.30, flap_reference_position=(1.10, -0.05),
).configure(flap_deflection_deg=25.0)
conditions = OperatingConditions(alpha_deg=6.0, reynolds=200000.0, mach=0.05)
runner = SolverRunner(Path.cwd(), "outputs/my_solver_attempts")
result = runner.run(MockSolver(), geometry, conditions,
                    min_clearance=0.002, timeout_s=0.5)
restored = SolverResult.load(runner.execution_root / result.attempt_directory / "result.json")
assert restored.to_dict() == result.to_dict()
assert not result.physical_coefficients_available
```

Uma nova chamada sempre cria outra tentativa. Não há cache de resultados,
paralelismo de análises, distribuição, correção automática de geometria nem
invocação de solver físico. Os números acima são entradas explícitas para
demonstrar infraestrutura, sem afirmar viabilidade aerodinâmica.

## Módulos e API pública

| Módulo | Responsabilidade |
|---|---|
| `conditions.py` | `OperatingConditions`, `ConditionParameter` |
| `base.py` | `SolverIdentity`, `SolverCase`, contrato `AerodynamicSolver` |
| `results.py` | Estados, coeficientes, referência de momento, diagnósticos, `SolverResult` |
| `execution.py` | `AttemptContext`, `ExecutionEvidence`, processo externo controlado |
| `runner.py` | Validação, etapas, isolamento, inventário e persistência |
| `serialization.py` | JSON finito, `FrozenJSON`, codificação canônica e SHA-256 |
| `mock.py` | Cenários determinísticos não físicos |
| `_paths.py` | Contenção de caminhos e compatibilidade de nomes Windows |
| `_windows_job.py`, `_process_worker.py` | Árvore de processos no Windows e porta de liberação |

Os tipos públicos são reexportados em `aeroopt.solvers`. Adaptadores não devem
usar os módulos privados; usem os métodos de `AttemptContext` para arquivos e
execução. Os contratos de dados são dataclasses imutáveis; arrays são snapshots
independentes em buffers sem escrita. `FrozenJSON` guarda uma representação
imutável e retorna cópias novas por `to_dict()`, inclusive para estruturas aninhadas.

## Condições, eixos e unidades

`OperatingConditions(alpha_deg, reynolds, mach, additional=())` exige os três
valores, sem padrões de escoamento. Reais finitos são aceitos; bools, strings,
NaN e infinitos são rejeitados. Reynolds é estritamente positivo e baseado na
**corda principal**; Mach é não negativo. Não há limites de alpha/Re/Mach
específicos de solver no núcleo.

A convenção explícita de alpha é a direção da velocidade a jusante no referencial
global: `V_inf / |V_inf| = (cos(alpha), sin(alpha))`. Assim alpha positivo tem
componente de velocidade positiva em y. O JSON registra `alpha_convention` e
`reynolds_reference=main_chord`. Não se gira a geometria para aplicar alpha. Um
adaptador deve converter seus sinais/eixos de forma documentada e verificável.
As posições globais continuam em `x/c_main,y/c_main`, com x a jusante, y para
cima e deflexão positiva do flap horária. Não se confunde o sinal de deflexão
com uma convenção de momento.

Parâmetros extras são uma tupla de `ConditionParameter(name, value, unit)`;
valores permitidos: float/int/string/bool. Extensões numéricas precisam de
unidade explícita, incluindo `dimensionless` quando apropriado. Nomes devem ser
únicos e são ordenados para serialização. O núcleo não inventa transição,
turbulência, densidade, viscosidade nem condições de contorno. O adaptador valida
quais extensões entende, seus domínios e a interpretação por elemento/superfície.

`AerodynamicCoefficients` admite `cl`, `cd`, `cm` ausentes (`None`/`null`). Zero é
um valor presente e não substitui ausência. Valores presentes precisam ser reais
finitos. Um adaptador real deve declarar e verificar eixos, sinais e normalização
das forças na corda principal, sem converter convenções silenciosamente.

`cm` presente exige `MomentReference(x_over_c, y_over_c, positive_sense,
normalization)`. A posição é global, em cordas principais; o sentido precisa
ser `clockwise` ou `counterclockwise`, e a definição de normalização é texto
obrigatório do adaptador. Não há referência implícita em 1/4 de corda, sinal
default ou transferência automática de momento. O núcleo não interpreta a
fórmula textual como uma conversão física.

## Etapas do adaptador

`AerodynamicSolver` define:

1. `identity`: identificação/versão, flag explícita de simulação, necessidade de
   convergência, capacidades e configurações `FrozenJSON`.
2. `validate_input(geometry, conditions) -> tuple[str, ...]`: rejeições específicas
   de faixa/modelo, após a validação geométrica comum. Tupla vazia significa aceitação.
3. `prepare(case, context)`: arquivos próprios do adaptador.
4. `execute(case, context) -> ExecutionEvidence`: execução da biblioteca ou
   `context.run_process(argv, ...)`.
5. `parse(case, context, evidence) -> ParsedOutput`: leitura dos resultados e
   declaração explícita de convergência/coeficientes.

O runner faz a persistência e captura falhas por etapa. `SolverCase` recebe a
montagem `MultiElementAirfoil` e as condições, mas fornece `main_points` e
`flap_points` globais congelados: estes são os dados a preparar/executar.
O nome/caminho da fonte e os parâmetros geométricos ficam na proveniência.
O runner verifica que manifesto e DAT globais continuam iguais após preparação
e após execução; arquivos específicos adicionais podem ter seus próprios formatos.
A identidade do adaptador também não pode mudar durante a tentativa.

## Estados e convergência

| `RunStatus` | Significado |
|---|---|
| `success` | Etapas completas e política de convergência atendida |
| `invalid_input` | Contrato, geometria/folga ou validação do adaptador rejeitou a entrada |
| `preparation_failed` | Falha ao preparar arquivos ou violação do manifesto durante preparação |
| `execution_failed` | Exceção de execução, saída não zero ou falha de integridade dos artefatos |
| `timeout` | Deadline externo expirado, com tentativa de encerramento da árvore |
| `parsing_failed` | Saída inválida, valores não finitos, dados globais alterados ou violação do contrato de parsing |
| `not_converged` | Convergência negativa ou não estabelecida segundo a política do adaptador |

`ConvergenceStatus` é separado: `unknown`, `not_applicable`, `converged` ou
`not_converged`. Código de saída zero não estabelece convergência. Um adaptador
com `convergence_required=True` precisa retornar `converged` para sucesso.
`not_applicable` pode produzir sucesso somente quando a identidade declara que
convergência não é exigida, como no mock. `unknown` nunca produz sucesso.

Coeficientes eventualmente lidos de uma análise real não convergida podem ser
preservados com o status correspondente, mas não ficam elegíveis como resultado
físico disponível. `physical_coefficients_available` requer sucesso, adaptador
declarado real e ao menos um coeficiente presente. O rótulo é derivado e
revalidado na reimportação, não um campo livre confiado ao JSON.

Para adaptadores simulados, qualquer `cl`, `cd` ou `cm` presente é proibido
mesmo que o número seja zero. A restrição vale no runner e em `SolverResult`,
inclusive ao carregar arquivos. O mock mantém os três `null` em todos os cenários.

## Identidade científica e tentativas

`SolverCase.case_id` é SHA-256 sobre uma árvore canônica contendo:

- Coordenadas float64 globais completas e separadas, na ordem das amostras.
- Referencial global e convenção de fechamento da fuga.
- Condições e extensões, incluindo referência de Reynolds e convenção de alpha.
- Identificação/versão/configurações/capacidades do adaptador e flag de simulação.
- Versão `aeroopt.solver-input.v1` e folga mínima exigida.

Antes de SHA-256, cada nó tem tag de tipo: objeto, array, float64, int, str,
bool ou null. Objetos são ordenados por chave; arrays preservam ordem. Float64
usa `float.hex()` exato; -0.0 e +0.0 são equivalentes apenas no hash. As tags
dos containers impedem que um dicionário fornecido imite um nó float. O texto
canônico é codificado em UTF-8. Isso não é uma implementação de RFC 8785.

Não há arredondamento/tolerância no hash: uma diferença de um bit pode mudar o
ID, assim como mudar a discretização ou a ordem dos pontos. Entradas numericamente
equivalentes segundo alguma interpretação externa não são automaticamente
agrupadas. Campos operacionais principais são convertidos explicitamente para
float64; tipos int/float das extensões/configurações permanecem distintos.
Texto Unicode não recebe normalização automática.

Caminhos de entrada/saída, nomes de exibição, timestamps, PID, duração e UUID
não entram na identidade. Timeout é política da tentativa, não da análise;
uma repetição com outro prazo conserva o ID. O adaptador precisa versionar sua
implementação e declarar configurações/bibliotecas/binário/ambiente que alterem
o significado científico. Dependências e variáveis de ambiente herdadas não
são automaticamente incorporadas ao ID: as versões de runtime registradas são
evidência, não garantia de identidade bit a bit entre plataformas.

Cada chamada ganha `attempt_id` UUID4 independente. Não há reutilização ou
sobrescrita de tentativa. Entradas impossíveis de representar canonicamente
(por exemplo transformação com overflow) ficam com `case_id=null`, sem inventar
um identificador científico, mas mantêm tentativa e diagnóstico próprios.

## Diretórios, arquivos e caminhos

```text
outputs/solver_interface_v020_final/
  summary.json                         # demonstração
  attempt-<UUID4>/
    attempt.json                       # alocação e política inicial
    input.json                         # snapshot científico
    geometry_validation.json
    geometry/main.dat
    geometry/flap.dat
    events.jsonl                       # histórico incremental de etapas
    mock_request.json
    adapter.log
    mock_output.json                   # também preservado quando inválido
    partial_prepare.txt                # cenário de falha de preparação
    traceback.txt                      # quando há exceção
    process.stdout.log                 # quando há processo externo
    process.stderr.log
    process.execution.json
    execution.json
    result.json                        # registro final e inventário
```

Só existem arquivos das etapas alcançadas. A estrutura de tentativas é **plana**:
o hash científico fica nos manifestos, não acrescenta 64 caracteres ao caminho.
Isso reduz problemas de executáveis Windows que usam APIs legadas com arquivos
relativos. O framework usa internamente caminhos longos `\\?\` no Windows;
programas externos podem continuar tendo seus próprios limites. Para eles,
escolha uma raiz curta e valide seu comportamento. Nenhuma junction/drive
virtual é criada para contornar limitações de executáveis.

`SolverRunner(workspace, directory)` exige uma raiz existente e um diretório
relativo nela. `context.path`, `write_text`, `write_json` e os rótulos de logs
rejeitam caminhos absolutos, `..`, nomes de dispositivos, alternate data streams,
symlinks e junctions. Logs/artefatos no resultado são relativos à tentativa;
`attempt_directory` é relativo à raiz de execução. Caminhos absolutos podem
aparecer apenas como evidência de origem/comando/runtime, fora da identidade.

Os adaptadores e executáveis são código **confiável**: essa contenção de caminhos
da API não é um sandbox de sistema operacional. Um binário malicioso poderia
ignorar a API e acessar outros arquivos. Não há proteção contra alterações
concorrentes hostis do filesystem ou autenticação criptográfica de resultados.

## Execução externa e timeout

`context.run_process([absolute_executable, arg1, ...], label="process")` recebe
uma sequência de argumentos e usa `shell=False`. No Windows exige `.exe`
explícito, sem execução implícita de `.bat/.cmd/.ps1`. Arquivos Python devem
usar explicitamente o interpretador. O diretório corrente é a tentativa;
TEMP/TMP/TMPDIR e o cache Matplotlib do processo apontam para seu `tmp/`.
O restante do ambiente é herdado e não é despejado em logs.

Stdout e stderr vão para arquivos, evitando bloqueio de pipes cheios.
`ExecutionEvidence` registra retorno, tempo, deadline expirado, comando,
caminhos de logs, PID do supervisor e evidência de encerramento da árvore.
Para bibliotecas Python, o código de processo pode ser ausente (`None`).

No Windows 8+/11, um worker Python isolado (`-I`) espera uma porta de liberação
em stdin. Antes de liberar o executável, o runner associa esse worker a um
Job Object com `JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE`, sem permissão de breakaway.
Assim, o executável não pode criar filhos na janela anterior à associação.
Se a associação falhar, o worker é encerrado sem liberar o programa-alvo.
Após saída normal, falha ou timeout, o job é terminado e seu contador de
processos ativos precisa chegar a zero; o supervisor é reaproveitado pelo SO
somente após `wait`. Não há uso de APIs privadas do CPython ou janelas visíveis.

Em POSIX há implementação por grupo de processos, mas esta entrega foi
executada em Windows, sem alegar validação equivalente nesse sistema. Filhos
comuns são abrangidos; mecanismos deliberados de escape de grupo/job ou criação
via serviços/WMI ficam fora do contrato de executáveis confiáveis.

O timeout opcional precisa ser positivo e finito. O deadline começa na etapa
de execução e restringe chamadas externas subsequentes; um adaptador pode pedir
um prazo menor, nunca ampliar o restante disponível. O prazo não interrompe
arbitrariamente preparação, parsing, chamadas bloqueantes em bibliotecas Python
ou a criação inicial de processo no SO. Limpeza tem limite próprio de 5 s.
Um adaptador que precise de deadline rígido deve executar a biblioteca em
processo externo controlado. `None` significa ausência de prazo configurado.

## Persistência, JSON e auditoria

`aeroopt.solver-result.v1` usa os seguintes grupos obrigatórios:

| Campos | Conteúdo |
|---|---|
| `case_id`, `attempt_id`, `attempt_directory` | Identidades separadas e caminho relativo |
| `solver`, `conditions` | Contratos completos; podem ser ausentes se nem foram construídos |
| `status`, `convergence` | Estados independentes |
| `coefficients`, `physical_coefficients_available` | Valores opcionais, referência de momento e elegibilidade derivada |
| `diagnostics` | Lista com `stage`, `code`, `message` |
| `duration_s`, `execution_duration_s` | Segundos monotônicos; execução ausente quando não alcançada |
| `started_at_utc`, `finished_at_utc` | Timestamps UTC conscientes de timezone |
| `artifacts` | Caminho relativo, SHA-256 dos bytes e tamanho |
| `provenance` | Versão AeroOpt, Python/dependências/plataforma, política e fontes geométricas |
| `input_payload`, `geometry_validation` | Snapshot que reproduz o hash e relatório geométrico estruturado |

A proveniência guarda hash dos valores numéricos originais imutáveis, índice/método
LE e parâmetros. O caminho da fonte é apenas de exibição: não se afirma que um
arquivo atualmente nesse caminho ainda tenha os bytes usados na importação.
O snapshot global e os DAT da tentativa são a entrada efetivamente auditável.

O JSON de evidências usa os decimais round-trip do Python, ordenação de chaves
e UTF-8. Não permite NaN/Infinity, números que transbordam em leitura ou chaves
duplicadas. A reimportação revalida tipos, estados, identidade científica,
simulação e a elegibilidade física. Os dataclasses tipados são o esquema
executável; não há dependência de um validador JSON Schema externo.

Publicação JSON escreve e sincroniza um arquivo temporário completo, depois
cria um hard link exclusivo para o destino e remove o temporário. Requer
filesystem com hard links (NTFS local verificado). Um destino existente nunca
é substituído. A publicação de cada JSON é atômica; a tentativa inteira não
é uma transação. Falta de espaço, perda de energia ou permissão que impeça
persistência pode deixar evidência parcial; `OSError` é propagado se não for
possível salvar o relatório, sem fingir sucesso de auditoria.

O inventário inclui arquivos parciais/saídas malformadas e logs, com hashes de
bytes. Exclui `result.json` para evitar autorreferência de hash. `events.jsonl`
é incremental, não uma transação. Exceções comuns e `SystemExit` de biblioteca
são convertidos em diagnóstico da etapa; `KeyboardInterrupt` é persistido e
depois relançado para respeitar a interrupção do operador. O encerramento do
processo externo ocorre em `finally`.

## Requisitos para um adaptador real

- Identidade/versão autênticas, capacidades, faixas suportadas, configurações e
  dependências cientificamente relevantes declaradas.
- Referências LE/TE revisadas, qualidade/discretização adequada ao solver e
  conversões de coordenadas/eixos verificadas, sem reparar entradas silenciosamente.
- Condições explícitas, incluindo transição/modelagem quando exigidas; rejeição
  explicativa de condições não suportadas.
- Critérios de convergência realmente lidos da saída, detecção de truncamentos,
  falhas numéricas e validade de coeficientes.
- Convenções de força e momento verificadas e conjuntos de referência
  independentes, incluindo casos de falha e sensibilidade à discretização.
- Uso do contexto para timeout/artefatos e testes do executável no Windows.

Uma saída de processo zero e geometria válida são pré-requisitos operacionais,
não validação física do modelo. Comparação entre fidelidades e integração à
otimização permanecem para etapas posteriores.

## Referências das APIs de execução

- [Python 3.13 — subprocess](https://docs.python.org/3.13/library/subprocess.html).
- [Microsoft — Job Objects](https://learn.microsoft.com/en-us/windows/win32/procthread/job-objects).
- [Microsoft — estruturas de limites de Job Object](https://learn.microsoft.com/en-us/windows/win32/api/winnt/ns-winnt-jobobject_extended_limit_information).

As definições de identidade e política científica acima são decisões do AeroOpt.
