# AI Code Agent rules

## Git Commit Convention

- For co-author, please add AI agent and model name as co-author.
  - Multi-co-author format (note no extra newline between co-authors):
    ```
    Co-authored-by: Claude Code <noreply@anthropic.com>
    Co-authored-by: glm-5 <service@zhipuai.cn>
    ```
  - First co-author: Agent
    - Claude Code: noreply@anthropic.com
  - Second co-author: Model
    - qwen* (eg. qwen3.5-plus): qianwen_opensource@alibabacloud.com
    - glm* (eg. glm-5): service@zhipuai.cn
    - minimax* (eg. MiniMax-M2.5): model@minimax.io
    - deepseek* (eg. DeepSeek-V3.2): service@deepseek.com
    - kimi* (eg. kimi-k2.5): growth@moonshot.cn
  - Model name should include the version or details, such as `qwen3.5-plus`, `glm-5`, which can be inferred by Claude Code's `/model` property.

