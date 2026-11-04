# Instruções permanentes — AeroOpt

- Consulte `docs/PROJECT_CONTEXT.md` antes de alterar a arquitetura.
- Preserve e consulte `docs/GEOMETRY_CONVENTIONS.md` para alterações geométricas.
- Use exclusivamente `.venv\Scripts\python.exe` para operações Python do projeto.
  Não recrie a `.venv` nem instale dependências globalmente.
- Execute os testes após modificações: `.\.venv\Scripts\python.exe -m pytest -q`.
- Atualize a documentação quando decisões técnicas ou a API forem alteradas.
- Não confunda resultados geométricos com resultados aerodinâmicos.
- Não modifique dados de entrada silenciosamente. Não suavize, interpole, reordene
  ou redistribua pontos automaticamente.
- Não introduza funcionalidades fora do escopo sem justificativa documentada.
- Não execute comandos destrutivos sem autorização. Preserve o histórico Git,
  não altere remotos e não crie commits automaticamente.
- Mantenha arquivos gerados, caches e temporários dentro do workspace.
- A v0.1 é exclusivamente geométrica. CFD, solvers aerodinâmicos, NSGA-II e
  modelos substitutos pertencem a versões futuras.
- A v0.1.1 é uma etapa de validação geométrica. Consulte também
  `docs/VALIDATION_V011.md` antes de alterar referências ou tolerâncias.
- Preserve os downloads originais com SHA-256 e proveniência. Não versione
  coordenadas externas sem condições de redistribuição estabelecidas; os
  originais UIUC locais e seus derivados permanecem em `outputs/`.
- Para incluir os testes publicados: `.\.venv\Scripts\python.exe -m pytest -q --reference-data=outputs/reference_data_v011`.
  Sem os dados, as verificações analíticas permanecem executáveis offline;
  nunca invente medições externas.
- Não confunda origem nominal NACA, mínimo x e referência por maior distância
  à fuga. Use índices explícitos quando houver referência conhecida.
- Contato/folga são globais; duplicação e degeneração são verificações
  intrínsecas na escala própria do elemento. Não relaxe tolerâncias para aprovar
  testes sem reproduzir e justificar o defeito.

